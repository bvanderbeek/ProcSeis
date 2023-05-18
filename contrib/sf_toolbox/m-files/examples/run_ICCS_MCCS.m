%% Automated picking of relative delay times
clear;
close all;
clc;

% Make faster by setting maximum lag and or trimming traces before
% cross-correlation and parallelising

% Requires LOMAX directory
addpath /Users/bvanderbeek/research/software/sf_toolbox/m-files
addpath /Users/bvanderbeek/research/software/sf_toolbox/m-files/utils
setenv('LOMAX','/Users/bvanderbeek/research/software/matlab_toolboxes/FilterPicker_Lomax/');

% Seismogram parameters
prjDir     = '/Users/bvanderbeek/research/NEWTON/AxiSEM_SPECFEM3D_Projects/projects/ISO_Mirrored';
channel    = {'BXX','BXY'}; % Channels to load
ext        = 'semv'; % File extension
tf_bin     = false; % True if file format is binary, false for ascii
dataWin    = [-100 100]; % Data window (s)
ttap       = 15; % Length of taper applied to trace ends...should be <= fmax (s)
filtopts   = [1/40 1/15 2 1]; % Butterworth filter parameters [f_min f_max order tf_zero_phase]
tf_rotate  = true; % Rotate seismograms into transverse, radial, and vertical components?
tf_moveout = true; % Subtract 1D travel-time from traces?
aPhase     = 'S'; % Phase on which seismograms are flattened/ Phase that will be picked
aModel   = 'iasp91'; % 1D Model name

% What channel are we picking?
chan_ind    = 1; %[1,2];

% Lomax Filter Picker Parameters
PickParams.method = 'lomax';
PickParams.stw    = 30;
PickParams.ltw    = 90;
PickParams.ptw    = 10;
PickParams.E1     = 10;
PickParams.E2     = 10;
PickParams.tadd   = 15;

% ICCS Parameters
iwstart = -5; % Start of cross-correlation window (s)
iwend   = 20; % End of cross-correlation window (s)
itaper  = 15; % Length of taper on window (s)
ccf_min = 0.7; % Reject any arrivals where the cross-correlation coefficient is less than this value

% MCC parameters
mwstart = -5; % Start of cross-correlation window (s)
mwend   = 15; % End of cross-correlation window (s)
mtaper  = 5;  % Length of taper on window (s)

%% Analyze Data

% Load stations
station = sf_read_stations(prjDir);

% Load events
event = sf_read_events(prjDir);

% Loop over events
for ievt = 1:length(event.id)
    % Load Seismograms
    if tf_bin
        data = sf_load_bin_seis(station,event,event.id(ievt),channel,prjDir,ext,dataWin,...
            'filtopts',filtopts,'tf_rotate',tf_rotate,'tf_moveout',tf_moveout,'aPhase',aPhase,'aModel',aModel,'ttap',ttap);
    else
        data = sf_load_ascii_seis(station,event,event.id(ievt),channel,prjDir,ext,dataWin,...
            'filtopts',filtopts,'tf_rotate',tf_rotate,'tf_moveout',tf_moveout,'aPhase',aPhase,'aModel',aModel,'ttap',ttap);
    end
    % Loop over channels
    if ~isempty(data.seis)
        tic
        for jchn = 1:length(chan_ind)
            % Select channel index
            ichan = chan_ind(jchn);
            
            % Iteratively align traces (ICCS)
            [delays,tpick1,ibad,CCF,itol,S1] = ICCS(data.t(1),data.fs,squeeze(data.seis(:,ichan,:)),itaper,iwstart,iwend-iwstart,ccf_min,PickParams);
            tpick1 = tpick1 + PickParams.tadd; % A user-prescribed ad-hoc adjustment
            % COULD ITERATE ON ICCS HERE
            % Clean bad data, re-stack and re-pick, adjust correlation window, run
            % ICCS alignment, repeat as needed
            
            % Get relative delay times via MCCC
            % Only running with good traces
            [DT,S2,ddt1] = MCCS(data.t(1),data.fs,squeeze(data.seis(~ibad,ichan,:)),tpick1+delays(~ibad),mtaper,mwstart,mwend-mwstart);
            
            % Final Pick on MCC stack
            % Lomax Filter Picker Method
            % Reduce tolerance until pick is made
            tpick2 = [];
            dE    = 0;
            E1    = PickParams.E1;
            E2    = PickParams.E2;
            while isempty(tpick2) && (E1 > 0) && (E2 > 0)
                % Filter picker parameters
                stw   = round(data.fs*PickParams.stw);
                ltw   = round(data.fs*PickParams.ltw);
                ptw   = round(data.fs*PickParams.ptw);
                E1    = PickParams.E1-dE;
                E2    = PickParams.E2-dE;
                % Make pick
                tpick2 = filterPicker(1/data.fs,S2(:)',[stw ltw ptw E1 E2],[1 data.nsamp]);
                tpick2 = data.t(1) + tpick2;
                % Update pick tolerances
                dE    = dE + 1;
            end
            if (dE - 1) > 0
                warning(['Tolerance reduced by ',num2str(dE-1),' before a pick was identified.'])
            end
            if isempty(tpick2)
                warning('No pick identified. Setting arrival to 0 s.');
                tpick2  = 0;
            elseif length(tpick2) > 1
                % This is an ad-hoc selection in the case that multiple picks are
                % found
                warning('Multiple phases identified. Chosing earliest arrival...')
                tpick2  = tpick2(1);
            end
            tpick2 = tpick2 + PickParams.tadd; % A user-prescribed ad-hoc adjustment
            
            
            % Make an error estimate based on similarity to stacked
            % waveform (e.g. Chevrot, 2002).
            w      = get_window(data.nsamp,data.fs,mtaper,data.t(1),tpick2+mwstart,mwend-mwstart);
            sref   = w(:).*S2(:); % Windowed reference waveform
            sref   = sref./rms(sref(logical(w))); % Amplitude normalized reference waveform
            AC     = xcorr(sref,sref);
            ac_max = max(AC);
            sid    = station.id(~ibad);
            CCF    = zeros(size(DT));
            ddt2    = zeros(size(DT));
            for nsta = 1:length(sid)
                ist      = strcmp(data.station,sid{nsta});
                w        = get_window(data.nsamp,data.fs,mtaper,data.t(1),tpick2+DT(nsta)+mwstart,mwend-mwstart);
                si       = w(:).*squeeze(data.seis(ist,ichan,:));
                si       = si./rms(si(logical(w)));
                [r,lags] = xcorr(sref,si);
                CCF(nsta) = min(max(r),ac_max);
                iddt      = find((AC >= CCF(nsta)),1,'last');
                % Chose nearest value or interpolate
                % sdt(nsta) = max(lags(AC >= CCF(nsta)))/data.fs;
                ddt2(nsta) = (((lags(iddt+1)-lags(iddt))/(AC(iddt+1)-AC(iddt)))*(CCF(nsta)-AC(iddt)) + lags(iddt))/data.fs;
            end
            
            % Save results
            picks.eventid    = event.id(ievt);
            picks.phase      = aPhase;
            picks.channel    = data.channel(ichan);
            picks.TRZ        = tf_rotate;
            picks.station    = station.id(~ibad);
            picks.traveltime = tpick2 + data.tt1D(~ibad) + DT;
            picks.dt1        = ddt1;
            picks.ddt2       = ddt2;
            picks.CCF        = CCF;
            picks.filtopts   = filtopts;
            % Store ICCS parameters
            picks.iccs.iwstart  = iwstart;
            picks.iccs.iwend    = iwend;
            picks.iccs.itaper   = itaper;
            picks.iccs.ccf_min  = ccf_min;
            % Store MCCS parameters
            picks.mccs.mwstart = mwstart;
            picks.mccs.mwend   = mwend;
            picks.mccs.mtaper  = mtaper;
            % Store picking parameters
            picks.PickParams = PickParams;
            save([prjDir,'/SEIS/CMTSOLUTION_',num2str(event.id(ievt)),'/picks_',aPhase,'_',data.channel{ichan},'.mat'],'picks');
            clear picks
        end
        toc
        disp(['Finished Event ',num2str(ievt),' of ',num2str(length(event.id)),'...']);
    else
        disp(['No data for Event ID ',num2str(event.id(ievt)),' Skipping...']);
    end
end
