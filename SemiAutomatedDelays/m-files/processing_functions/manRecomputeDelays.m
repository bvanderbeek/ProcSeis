function [ATT,data] = manRecomputeDelays(evtid)

% Define analysis parameters
tf_specfem   = false; % (true) SPECFEM binary seismograms; (false) mseed seismograms
dataDir      = '/Users/bvanderbeek/research/CASCADIA/Waveforms_S'; % Data directory
theEvent     = [dataDir,'/event_cascadia.mat']; % Event structure
theStation   = [dataDir,'/station_cascadia.mat']; % Station structure
theATT       = [dataDir,'/ATT_cascadia.mat']; % The arrival time table
theBound     = [dataDir,'/PB_Cascadia.mat']; % Plate boundaries structure
ichan        = 1; % Initial channel to pick (1 = T, 2 = R,Q, 3 = Z,L)
sampFreq     = 40; % Desired Sampling frequency (Hz)
timeWindow   = 600; % Length of seismograms to load (s)
refModel     = 'ak135'; % Reference model for 1D travel-time predictions
SeismicPhase = 'S'; % A seismic phase to analyse
filtType     = 0; % Filter type (0 = Butterworth; only option currently implemented)
corners      = [1/33, 1/12]; % Corner frequencies for Butterworth filter (Hz)
order        = 3; % Order of Butterworth filter
tf_zerophase = true; % Use zero-phase filter?
% If true, initial delays are set to existing picks. If false, initial
% delays are assumed to be zero.
tf_initialAlignment = false;
% Analysis window parameters. Two time windows.
twin1min = -5;
twin1max = 15;
twin2min = 0;
twin2max = 10;
% Polarisation analysis window
pmin     = -5;
pmax     = 20;
tf_pelv = true; % Include elevation when rotating into PAZ?
% Multi-channel Cross-correlation parameters
nmax_mcc        = 20; % Maximum number of iterations
tol_mcc         = 2/sampFreq; % Tolerance for MCC (exit when delays change by less than tol)
tf_weighted_mcc = false; % Use correlation-coefficient weighted MCC?
% Auto save updated arrival time table at end?
tf_autoSave = true;
% Set environment variables
setenv('TAUPJAR','/Users/bvanderbeek/research/software/TauP_Toolkit/TauP-2.4.5/lib/TauP-2.4.5.jar');
% Add required paths
addpath /Users/bvanderbeek/research/SemiAutomatedDelays/contrib/mk_rd_mseed/
% Figure directories
if tf_specfem
    tmpDir   = [dataDir,'/SEIS/CMTSOLUTION_',num2str(evtid),'/tmpFigs'];
    figDir   = [dataDir,'/SEIS/CMTSOLUTION_',num2str(evtid)];
else
    tmpDir   = [dataDir,'/EVT',num2str(evtid),'/tmpFigs'];
    figDir   = [dataDir,'/EVT',num2str(evtid)];
end
% Make temporary figure directory
[status,msg] = mkdir(tmpDir); %#ok No warning for existing directory

% Load data structures
load(theEvent,'event');
load(theStation,'station');
load(theATT,'ATT');
if isempty(theBound)
    PB.boundaries = zeros(0,2);
else
    load(theBound,'PB');
end

% Get phase
aPhase = unique(ATT.phase(ATT.event == evtid));
if length(aPhase) == 1
    aPhase = aPhase{1};
elseif isempty(aPhase)
    aPhase = SeismicPhase;
else
    if evtid == 697
        warning('Assuming S for event 697.');
        aPhase = 'S';
    else
        error('Multiple phases recorded by this event!');
    end
end

% Load filtered seismograms
if tf_specfem
    sfEvent = sf_read_events(dataDir);
    sfStation = sf_read_stations(dataDir);
    data = sf_load_bin_seis(sfStation,sfEvent,evtid,{'BXX','BXY','BXZ'},dataDir,'semv',timeWindow*[-0.5, 0.5],...
        'filtopts',[corners(1),corners(2),order,tf_zerophase],'tf_rotate',true,...
        'tf_moveout',true,'ttap',[],'aPhase',aPhase,'aModel',refModel);
    % Modify data structure
    data.t = data.t(:)';
    data = rmfield(data,'event');
    data.event.id         = evtid;
    data.event.originTime = datenum(event.originDate(event.id == evtid)); % DATENUM
    data.event.latitude   = event.latitude(event.id == evtid);
    data.event.longitude  = event.longitude(event.id == evtid);
    data.event.depth      = event.depth(event.id == evtid);
    if data.ntrace == length(station.name)
        data.nsta = data.ntrace;
        data = rmfield(data,'ntrace');
    else
        error('Missing traces!');
    end
    [~,lis] = ismember(data.station,station.name);
    data.network = station.network(lis);
    data.refTime = data.event.originTime + (data.tt1D./(60*60*24));
    data.seis = permute(data.seis,[1,3,2]);
else
    % Loads data into TRZ components
    data = loadSeis_commonEvent(dataDir,evtid,event,station,sampFreq,...
        timeWindow,refModel,aPhase,corners,order,tf_zerophase);
end

% Define analysis windows for ALIGNED traces
awin_1 = (data.t >= twin1min) & (data.t <= twin1max);
awin_2 = (data.t >= twin2min) & (data.t <= twin2max);
pwin   = (data.t >= pmin) & (data.t <= pmax);
% Define time limits for plots
tlim = [-max(1./corners), twin1max + max(1./corners)];

%% Associate Original Picks
figFile = [tmpDir,'/','0_PriorAlignment_Transverse'];

% Arrival identification parameters
CHN = 'T'; % Original channel picked
FLT = [0,1/33,1/12,3,0]; % Original filter used
WND = 0; % Only original picks should have a window parameter of 0 s

% Loop over stations and find picks
DT_i  = zeros(data.nsta,1);
ddt_i = zeros(data.nsta,1);
for ii = 1:data.nsta
    % Index arrival
    keep = (ATT.event == evtid) & strcmp(ATT.station,data.station{ii})...
        & strcmp(ATT.channel,CHN) & strcmp(ATT.phase,aPhase)...
        & ismember(ATT.filter,FLT,'rows') & (ATT.window == WND);
    
    % Evaluate matches
    if sum(keep) == 1
        AT        = datenum(ATT.arrivalTime(keep));
        DT_i(ii)  = (AT - data.refTime(ii))*24*60*60;
        ddt_i(ii) = ATT.error(keep);
    elseif sum(keep) == 0
        warning('No arrivals found.');
    else
        warning('Multiple arrivals found');
    end
end
meanDelay_i = mean(DT_i);
DT_i        = DT_i - meanDelay_i;

% Use original picks to initially align seismograms?
if tf_initialAlignment
    data.meanDelay = meanDelay_i;
    data.DT        = DT_i;
    data.ddt       = ddt_i;
else
    data.meanDelay = 0;
    data.DT        = zeros(data.nsta,1);
    data.ddt       = -ones(data.nsta,1);
end

% Align traces
S  = align_traces(data.t,squeeze(data.seis(:,:,ichan)),data.meanDelay + DT_i);
% Scale amplitudes
S  = scale_traces(awin_1,S,1);
S0 = align_traces(data.t,S,-(data.meanDelay + DT_i));
% Stack traces
D  = mean(S,1);
D0 = mean(S0,1);
% Plot
H = figure;
subplot(1,2,1); hold on;
plot(data.t,S0);
plot(data.t,D0,'-k','linewidth',2);
plot([twin1min*ones(2,1), twin1max*ones(2,1)],...
    max(abs(D(:)))*[-1,-1; 1,1],'-g','linewidth',2);
box on; grid on; axis tight; xlim(tlim);
subplot(1,2,2); hold on;
plot(data.t,S);
plot(data.t,D,'-k','linewidth',2);
plot([twin1min*ones(2,1), twin1max*ones(2,1)],...
    max(abs(D(:)))*[-1,-1; 1,1],'-g','linewidth',2);
box on; grid on; axis tight; xlim(tlim);
H.Position(3) = 2*H.Position(3);
% Pause for display
fprintf('\n Prior alignment. Press any key to continue. \n');
figure(H);
pause

% Option to Update Mean delays
r = 10;
while ~(r == 0) && ~(r == 1)
    r = input(' Update mean delay/shift window? (yes = 1; no = 0)? ');
end
if r == 0
    fprintf('\n Mean delay NOT updated. \n');
else
    fprintf('\n Updated mean delay. \n');
    % Pick stacked seismogram with new alignment
    data.meanDelay = pickStack(data.t,squeeze(data.seis(:,:,ichan)),data.meanDelay,data.DT,awin_1,tlim);
end

% Save initial alignment figure
figure(H);
print(figFile,'-dpng');
close(H);

% For later plotting
[~,ista]    = ismember(data.station,station.name);
priorDelays = [station.longitude(ista),station.latitude(ista),DT_i];

%% Check for Polarity Reversals

% Align traces
S = align_traces(data.t,squeeze(data.seis(:,:,ichan)),data.meanDelay + data.DT);
% Scale amplitudes
S = scale_traces(awin_1,S,1);
% Stack traces
D = mean(S,1);

% Zero-lag correlation with stacked trace
iflip = S(:,awin_1)*(D(:,awin_1)')./sqrt(sum(S(:,awin_1).^2,2).*sum(D(:,awin_1).^2,2));
% Allow for some anti-correlation
iflip = iflip < mean(iflip - 1); % Average distance from perfect correlation
% iflip = iflip < -std(iflip); % Standard-deviation in correlation
% Store flipped station indices for future plotting
stflip = ista(iflip);

% Display number of estimated polarity reversals
if sum(iflip) > 0
    fprintf(['Estimated ',num2str(sum(iflip)),' polarity reversals. \n']);
    fprintf('\n Stations: \n');
    disp(data.station(iflip)');
    
    % Plot
    H = figure; hold on;
    plot(data.t,S(~iflip,:));
    plot(data.t,D,'-k','linewidth',2);
    plot(data.t,S(iflip,:),'-r','linewidth',1);
    plot([twin1min*ones(2,1), twin1max*ones(2,1)],...
        max(abs(D(:)))*[-1,-1; 1,1],'-g','linewidth',2);
    box on; grid on; axis tight; xlim(tlim);
    
    % Option to flip polarities of suspected traces
    r = 10;
    while ~(r == 0) && ~(r == 1)
        r = input(' Change polarity of suspected traces (yes = 1; no = 0)? ');
    end
    
    % Evaluate response
    if r == 0
        fprintf('\n Polarities are NOT changed. \n');
    elseif r == 1
        fprintf('\n Selected trace polarities modified. Press any key to continue. \n');
        
        % Reverse polarity of selected trace
        data.seis(iflip,:,ichan) = -data.seis(iflip,:,ichan);
        
        % Re-stack
        S(iflip,:) = -S(iflip,:);
        D          = mean(S,1);
        
        % Re-plot
        close(H);
        H = figure; hold on;
        plot(data.t,S(~iflip,:));
        plot(data.t,D,'-k','linewidth',2);
        plot(data.t,S(iflip,:),'-r','linewidth',1);
        plot([twin1min*ones(2,1), twin1max*ones(2,1)],...
            max(abs(D(:)))*[-1,-1; 1,1],'-g','linewidth',2);
        box on; grid on; axis tight; xlim(tlim);
        pause;
    else
        error('Problem evaluating user response!');
    end
    
    % Close plot
    close(H);
end

%% 1) First Iterative Multi-channel Cross-correlation
figFile = [tmpDir,'/','1_FirstAlignment_Transverse'];

% Call iterative MCC
fprintf('\n Starting first iterative MCC...\n');
[DT,ddt,ibad,H] = call_iterateMCC(data.t,squeeze(data.seis(:,:,ichan)),data.fs,...
    data.meanDelay,data.DT,awin_1,twin1min,twin1max,...
    tf_weighted_mcc,tol_mcc,nmax_mcc,tlim);
fprintf(' ...Completed first iterative MCC \n');
fprintf([' Minimum/Maximum Relative Delays (ms): ',num2str(round(1000*min(DT))),'/',num2str(round(1000*max(DT))),'\n']);

% Evaluate alignment
[data,flag] = continueAlignment(data,DT,ddt,ibad,H,figFile);
if flag == 0
    fprintf('\n Exiting program. \n');
    return
end

% Pick stacked seismogram with new alignment
data.meanDelay = pickStack(data.t,squeeze(data.seis(:,:,ichan)),data.meanDelay,data.DT,awin_1,tlim);

%% 3) Second Iterative Multi-channel Cross-correlation
figFile = [tmpDir,'/','2_SecondAlignment_Transverse'];

if ~((twin2min == twin1min) && (twin2max == twin1max))
    % Call iterative MCC
    fprintf('\n Starting second iterative MCC...\n');
    [DT,ddt,ibad,H] = call_iterateMCC(data.t,squeeze(data.seis(:,:,ichan)),data.fs,...
        data.meanDelay,data.DT,awin_2,twin2min,twin2max,...
        tf_weighted_mcc,tol_mcc,nmax_mcc,tlim);
    fprintf(' ...Completed second iterative MCC \n');
    fprintf([' Minimum/Maximum Relative Delays (ms): ',num2str(round(1000*min(DT))),'/',num2str(round(1000*max(DT))),'\n']);
    
    % Evaluate alignment
    [data,flag] = continueAlignment(data,DT,ddt,ibad,H,figFile);
    if flag == 0
        fprintf('\n Exiting program. \n');
        return
    end
end

%% Update Arrival Time Table

CHN = data.channel; % {'T','R','Z'};
POL = {[0,0],[90,0],[0,90]};
FLT = [filtType,corners,order,~tf_zerophase];
ATT = update_ATT(ATT,data,CHN{ichan},aPhase,POL{ichan},FLT,twin2max-twin2min);

% Copy of delays and channel name for future plotting
[~,ista]  = ismember(data.station,station.name);
trzDelays = [station.longitude(ista),station.latitude(ista),data.DT];
if ichan == 1
    channelName = 'Transverse';
elseif ichan == 2
    channelName = 'Radial';
elseif ichan == 3
    channelName = 'Vertical';
else
    error('Channel index out-of-bounds!');
end

%% Check Polarity of Secondary Channels
% Just because the picked channel is polarity reversed doesn't imply the
% other channels are reversed. Now that we have a good initial alignment,
% let's check the polarity of the secondary channels before doing
% polarisation analysis.

if ~tf_specfem
    % Select non-picked channels
    jchan = [1,2,3];
    jchan(ichan) = [];
    
    % Align traces
    S1 = align_traces(data.t,squeeze(data.seis(:,:,jchan(1))),data.meanDelay + data.DT);
    S2 = align_traces(data.t,squeeze(data.seis(:,:,jchan(2))),data.meanDelay + data.DT);
    % Scale amplitudes
    S1 = scale_traces(awin_1,S1,1);
    S2 = scale_traces(awin_1,S2,1);
    % Stack traces
    D1 = mean(S1,1);
    D2 = mean(S2,1);
    % Zero-lag correlation with stacked trace
    iflip_1 = S1(:,awin_1)*(D1(:,awin_1)')./sqrt(sum(S1(:,awin_1).^2,2).*sum(D1(:,awin_1).^2,2));
    iflip_1 = iflip_1 < min(mean(iflip_1 - 1), -std(iflip_1));
    iflip_2 = S2(:,awin_1)*(D2(:,awin_1)')./sqrt(sum(S2(:,awin_1).^2,2).*sum(D2(:,awin_1).^2,2));
    iflip_2 = iflip_2 < min(mean(iflip_2 - 1), -std(iflip_2));
    % Plot check
    H = figure;
    subplot(2,1,1); hold on;
    plot(data.t,S1);
    plot(data.t,S1(iflip_1,:),'-r','linewidth',2);
    plot(data.t,D1,'-k','linewidth',2);
    plot([twin1min*ones(2,1), twin1max*ones(2,1)],...
        max(abs(D1(:)))*[-1,-1; 1,1],'-g','linewidth',2);
    box on; grid on; axis tight; title(num2str(jchan(1))); xlim(tlim);
    subplot(2,1,2); hold on;
    plot(data.t,S2);
    plot(data.t,S2(iflip_2,:),'-r','linewidth',2);
    plot(data.t,D2,'-k','linewidth',2);
    plot([twin1min*ones(2,1), twin1max*ones(2,1)],...
        max(abs(D2(:)))*[-1,-1; 1,1],'-g','linewidth',2);
    box on; grid on; axis tight; title(num2str(jchan(2))); xlim(tlim);
    figure(H);
    
    % Option to flip polarities of suspected traces
    r = 10;
    while ~(r == 0) && ~(r == 1)
        r = input(' Change polarity of suspected traces (yes = 1; no = 0)? ');
    end
    if r == 1
        data.seis(iflip_1,:,jchan(1)) = -data.seis(iflip_1,:,jchan(1));
        data.seis(iflip_2,:,jchan(2)) = -data.seis(iflip_2,:,jchan(2));
    end
    close(H);
end

%% Polarisation Analysis
figFile = [tmpDir,'/','3_PolarisationResults'];

% Align traces
TRZ = align_traces(data.t,data.seis,data.meanDelay + data.DT);
% Scale amplitudes
TRZ = scale_traces(pwin,TRZ,[]);
% Stack traces
TRZ = squeeze(mean(TRZ,1))';

% Polarisation analysis
[azm,elv] = compute_polarization(TRZ(1,pwin),TRZ(2,pwin),TRZ(3,pwin));
if tf_pelv
    PML = roty(elv)*rotz(-azm)*TRZ;
else
    elv = 0;
    PML = rotz(-azm)*TRZ;
end

% Evaluate polarisation alignment
H = figure;
% TRZ
subplot(1,2,1); hold on;
plot(data.t,TRZ(1,:),'-b','linewidth',2);
plot(data.t,TRZ(2,:),'-r','linewidth',2);
plot(data.t,TRZ(3,:),'-g','linewidth',2);
plot([pmin*ones(2,1), pmax*ones(2,1)],...
    max(abs(PML(:)))*[-1,-1; 1,1],'--k','linewidth',2);
box on; grid on; axis tight; xlim(tlim);
legend('transverse','radial','vertical');
title('TRZ');
% PML
subplot(1,2,2); hold on;
plot(data.t,PML(1,:),'-b','linewidth',2);
plot(data.t,PML(2,:),'-r','linewidth',2);
plot(data.t,PML(3,:),'-g','linewidth',2);
plot([pmin*ones(2,1), pmax*ones(2,1)],...
    max(abs(PML(:)))*[-1,-1; 1,1],'--k','linewidth',2);
box on; grid on; axis tight; xlim(tlim);
legend('polarisation','minor','logitudinal');
title(['Azimuth = ',num2str(round(azm)),', Elevation = ',num2str(round(elv))]);
% Double width
H.Position(3) = 2*H.Position(3);

% Evaluate results
r = 10;
while ~(r == 0) && ~(r == 1)
    r = input(' Continue with polarisation (yes = 1; no = 0)? ');
end
if r == 1
    % Update delays
    data.pol = [azm,elv];
    % Save figure
    figure(H);
    print(figFile,'-dpng');
    % Close figure and continue
    close(H);
elseif r == 0
    % Exit function
    fprintf('\n Exiting program. \n');
    return
else
    error('Problem evaluating user response!');
end

%% Polarisation: First Iterative Multi-channel Cross-correlation
figFile = [tmpDir,'/','4_FirstAlignment_Polarisation'];

% Define principle S-wave
Sp = zeros(data.nsta,data.nsamp);
for ii = 1:data.nsta
    % Rotate seismograms
    si = squeeze(data.seis(ii,:,:))';
    si = roty(elv)*rotz(-azm)*si;
    % Store principle component
    Sp(ii,:) = si(1,:);
end

% Pick stacked seismogram with new alignment
data.meanDelay = pickStack(data.t,Sp,data.meanDelay,data.DT,awin_1,tlim);

% Call iterative MCC
fprintf('\n Polarisation: starting first iterative MCC...\n');
[DT,ddt,ibad,H] = call_iterateMCC(data.t,Sp,data.fs,data.meanDelay,data.DT,...
    awin_1,twin1min,twin1max,tf_weighted_mcc,tol_mcc,nmax_mcc,tlim);
fprintf(' ...Completed first iterative MCC \n');
fprintf([' Minimum/Maximum Relative Delays (ms): ',num2str(round(1000*min(DT))),'/',num2str(round(1000*max(DT))),'\n']);

% Evaluate alignment
[data,flag] = continueAlignment(data,DT,ddt,ibad,H,figFile);
if flag == 0
    fprintf('\n Exiting program. \n');
    return
end

% Delete bad traces from principle coordinate seismograms
if sum(ibad) > 0
    Sp(ibad,:) = [];
end

% Pick stacked seismogram with new alignment
data.meanDelay = pickStack(data.t,Sp,data.meanDelay,data.DT,awin_1,tlim);

%% Polarisation: Second Iterative Multi-channel Cross-correlation
figFile = [tmpDir,'/','5_SecondAlignment_Polarisation'];

if ~((twin2min == twin1min) && (twin2max == twin1max))
    % Call iterative MCC
    fprintf('\n Polarisation: starting second iterative MCC...\n');
    [DT,ddt,ibad,H] = call_iterateMCC(data.t,Sp,data.fs,data.meanDelay,data.DT,...
        awin_2,twin2min,twin2max,tf_weighted_mcc,tol_mcc,nmax_mcc,tlim);
    fprintf(' ...Completed second iterative MCC \n');
    fprintf([' Minimum/Maximum Relative Delays (ms): ',num2str(round(1000*min(DT))),'/',num2str(round(1000*max(DT))),'\n']);
    
    % Evaluate alignment
    [data,flag] = continueAlignment(data,DT,ddt,ibad,H,figFile);
    if flag == 0
        fprintf('\n Exiting program. \n');
        return
    end
    
    % Delete bad traces from principle coordinate seismograms
    if sum(ibad) > 0
        Sp(ibad,:) = [];
    end
end

%% Plot Delays
figFile = [tmpDir,'/','6_RelativeDelays'];

% Station indexing
[~,ista] = ismember(data.station,station.name);
DT_max = max(abs(cat(1,priorDelays(:,3),trzDelays(:,3),data.DT(:))));

% Plot
H = figure;
% Prior Delays
subplot(1,3,1); hold on;
plot(PB.boundaries(:,1),PB.boundaries(:,2),'-k','linewidth',2);
scatter(priorDelays(:,1),priorDelays(:,2),75,priorDelays(:,3),'filled');
axis tight; axis image; box on; grid on; title('Initial');
colormap(jet); caxis(DT_max*[-1,1]); colorbar('southoutside');
% Transverse Delays
subplot(1,3,2); hold on;
plot(PB.boundaries(:,1),PB.boundaries(:,2),'-k','linewidth',2);
scatter(trzDelays(:,1),trzDelays(:,2),75,trzDelays(:,3),'filled');
axis tight; axis image; box on; grid on; title(channelName);
colormap(jet); caxis(DT_max*[-1,1]); colorbar('southoutside');
% Polarisation Delays
subplot(1,3,3); hold on;
plot(PB.boundaries(:,1),PB.boundaries(:,2),'-k','linewidth',2);
scatter(station.longitude(ista),station.latitude(ista),75,data.DT,'filled');
% Plot polarity flips
plot(station.longitude(stflip),station.latitude(stflip),'xk','markersize',10);
axis tight; axis image; box on; grid on; title('Polarisation');
colormap(jet); caxis(DT_max*[-1,1]); colorbar('southoutside');
% Increase figure width
H.Position(3) = 2*H.Position(3);

% Pause
fprintf('\n Relative delay comparison. Press any key to continue \n');
figure(H);
pause

% Save figure
figure(H);
print(figFile,'-dpng');

%% Update Arrival Time Table

CHN = 'PZ';
POL = data.pol;
FLT = [filtType,corners,order,~tf_zerophase];
ATT = update_ATT(ATT,data,CHN,aPhase,POL,FLT,twin2max-twin2min);

% Auto-save arrival time table
if tf_autoSave
    fprintf(['\n Saved arrival time table: ',theATT,'\n']);
    save(theATT,'ATT');
    
    % Move Figures to Data Directory
    figFiles = {'0_PriorAlignment_Transverse', '1_FirstAlignment_Transverse',...
        '2_SecondAlignment_Transverse', '3_PolarisationResults', '4_FirstAlignment_Polarisation',...
        '5_SecondAlignment_Polarisation', '6_RelativeDelays'};
    for ii = 1:length(figFiles)
        cmd = ['mv ', tmpDir, '/', figFiles{ii}, '.png ', figDir, '/.'];
        system(cmd);
    end
end

% Done!
fprintf(['\n Finished Event ',num2str(evtid),'! \n']);

end

%% subFunction: Update Arrival Time Table
function ATT = update_ATT(ATT,data,CHN,PHS,POL,FLT,WND)

% Create new table entries
EVT = data.event.id;
OT  = datetime(data.event.originTime,'ConvertFrom','datenum','format','yyyy-MM-dd HH:mm:ss.SSS');
NW  = data.network;
ST  = data.station;
% Arrival time (DATETIME)
AT  = data.refTime + ((data.meanDelay + data.DT)/(60*60*24));
AT  = datetime(AT,'ConvertFrom','datenum','format','yyyy-MM-dd HH:mm:ss.SSS');
ERR = data.ddt;
% Time of latest update
LUD = datetime;

% Loop over stations
for ii = 1:data.nsta
    
    % New table row
    ROW = {EVT,OT,NW{ii},ST{ii},CHN,PHS,POL,FLT,WND,AT(ii),ERR(ii),LUD};
    
    % Index arrival
    iup = (ATT.event == data.event.id) & strcmp(ATT.network,data.network{ii})...
        & strcmp(ATT.station,data.station{ii}) & strcmp(ATT.channel,CHN)...
        & strcmp(ATT.phase,PHS) & ismember(ATT.filter,FLT,'rows')...
        & (ATT.window == WND);
    
    % Evaluate matches
    if sum(iup) == 0
        % No match found; add new entry
        ATT = cat(1,ATT,ROW);
    elseif sum(iup) == 1
        % Unique match found; overwrite entry
        ATT(iup,:) = ROW;
    else
        % Multiple matches found; error
        error('Multiple matches!');
    end
end

end

%% subFunction: Call Iterative MCC
function [DT,ddt,ibad,H] = call_iterateMCC(t,S,fs,tphase,DT_i,w,tpre,tpost,...
    tf_weighted_mcc,tol_mcc,nmax_mcc,tlim)

% Call interative MCC
[DT,ddt] = iterateMCC(t,S,fs,tphase,DT_i,tpre,tpost,tf_weighted_mcc,tol_mcc,nmax_mcc);

% Check for anomalously large delays
ibad = abs(DT) > 8;

% Align traces
S0 = align_traces(t,S,tphase + DT_i);
S  = align_traces(t,S,tphase + DT);
% Scale amplitudes
S0 = scale_traces(w,S0,1);
S  = scale_traces(w,S,1);
% Stack traces
D0 = mean(S0,1);
D  = mean(S,1);

% Compare results
H = figure;
% Initial alignment
subplot(1,2,1); hold on;
plot(t,S0);
plot(t,D0,'-k','linewidth',2);
plot([tpre*ones(2,1), tpost*ones(2,1)],...
    max(abs(D0(:)))*[-1,-1; 1,1],'-g','linewidth',2);
box on; grid on; axis tight; xlim(tlim);
xlabel('reduced time (s)'); ylabel('amplitude'); title('Initial');
% New alignment
subplot(1,2,2); hold on;
plot(t,S);
plot(t,D,'-k','linewidth',2);
plot([tpre*ones(2,1), tpost*ones(2,1)],...
    max(abs(D(:)))*[-1,-1; 1,1],'-g','linewidth',2);
box on; grid on; axis tight; xlim(tlim);
xlabel('reduced time (s)'); ylabel('amplitude'); title('After MCC');
% Plot bad delays
if sum(ibad) > 0
    fprintf(['\n There are ',num2str(sum(ibad)),' delays > 8 s! These will be deleted. \n']);
    subplot(1,2,1);
    plot(t,S0(ibad,:),'-r','linewidth',1);
    subplot(1,2,2);
    plot(t,S(ibad,:),'-r','linewidth',1);
end
% Double the figure width
H.Position(3) = 2*H.Position(3);
end

%% subFunction: Continue with New Alignment
function [data,r] = continueAlignment(data,DT,ddt,ibad,H,figFile)

% Delete any bad picks
if sum(ibad) > 0
    DT(ibad)    = [];
    ddt(ibad,:) = [];
    
    data.network(ibad)  = [];
    data.station(ibad)  = [];
    data.nsta           = sum(~ibad);
    data.delta(ibad)    = [];
    data.baz(ibad)      = [];
    data.refTime(ibad)  = [];
    data.seis(ibad,:,:) = [];
    data.DT(ibad)       = [];
    data.ddt(ibad)      = [];
end


% Evaluate results
r = 10;
while ~(r == 0) && ~(r == 1)
    r = input(' Continue with new alignment (yes = 1; no = 0)? ');
end
if r == 1
    % Update delays
    data.DT  = DT;
    data.ddt = ddt(:,2);
    % Print figure
    figure(H);
    print(figFile,'-dpng');
    % Close figure and continue
    close(H);
    fprintf('\n Delays updated. \n');
elseif r == 0
    % Exit function
    return
else
    error('Problem evaluating user response!');
end
end

%% subFunction: Pick Stacked Trace
function tphase = pickStack(t,S,tphase,DT,w,tlim)
% Align traces
S = align_traces(t,S,tphase + DT);
% Scale amplitudes
S = scale_traces(w,S,1);
% Stack traces
D = mean(S,1);

fprintf('\n Pick new mean delay. \n');

% Plot stacked waveform
H = figure; hold on;
plot(t,D,'-k','linewidth',2);
box on; grid on; axis tight; xlim(tlim);

% Evaluate pick
r = 10;
while ~(r == 1)
    % Activate figure
    figure(H);
    
    % Make pick
    [dtp,~] = ginput(1);
    % Plot pick
    G = plot(dtp*ones(1,2),max(abs(D))*[-1,1],'-r','linewidth',2);
    
    % Keep pick?
    while ~(r == 0) && ~(r == 1)
        r = input(' Keep mean delay (yes = 1; no = 0)? ');
    end
    
    % Evaluate response
    if r == 0
        r = 10;
        delete(G);
    elseif r == 1
        % Update mean delay
        tphase = tphase + dtp;
    end
end
close(H);
end