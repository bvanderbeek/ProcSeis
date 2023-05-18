function [data,stlist] = sf_load_su_seis(station,event,evtid,channel,prjDir,ext,dataWin,varargin)

% WEIRD PROBLEM WITH SU FORMAT. A SEEMINGLY RANDOM DISTRIBUTION OF STATIONS
% ARE NOT SAVED

% Load the seismograms
% Create common time vector
% Resample trace to consistent frequency
% Interpolate trace to common time vector
% --> Ideally absolute or travel-time

% Options
% Sampling rate (default to that in first file)
% Trace length (default to that in first file)
% Taper width (default is 5%)
% If resampling is required, it cannot handle non-integer sampling frequencies

% NEED TO CHECK THAT TIMING DEFINITION IS UNIVERSAL

% Make sure we can find the toolbox for reading SU data format
if isempty(getenv('SEGYMAT'))
    error('Missing SEGYMAT environment variable that points to segyMAT toolbox for reading SU data.')
else
    addpath(getenv('SEGYMAT'));
end

% Define variable input arguments
var_name      = {'filtopts','tf_rotate','tf_moveout','ttap','aPhase','aModel'};
var_value     = {[],false,false,2,[],[]};
[~,ib]        = ismember(varargin(1:2:(length(varargin)-1)),var_name);
var_value(ib) = varargin(2:2:length(varargin)); %#ok used in following loop
for ivar = 1:length(var_name)
    eval([var_name{ivar},' = var_value{',num2str(ivar),'};']);
end

% Check optional input parameters
if tf_moveout && (isempty(aPhase) || isempty(aModel))
    error('Need to input aPhase and a 1D model to use moveout correction.');
end

% Get seismogram start times. This is defined differently for an AxiSEM
% coupled simulation.
tf_coupled = sf_read_input([prjDir,'/DATA/Par_file'],'COUPLE_WITH_INJECTION_TECHNIQUE','=',false);
ctype      = sf_read_input([prjDir,'/DATA/Par_file'],'INJECTION_TECHNIQUE_TYPE','=',true);
if strcmp(tf_coupled,'.true.') && (ctype == 2)
    tstart = dlmread([prjDir,'/','/source_files/reformat_CMTSOLUTION_',num2str(evtid),'.par']);
    tstart = tstart(2,1);
    % SPECFEM seismograms start at -2x the source duration.
    thd     = sf_read_input([prjDir,'/source_files/CMTSOLUTION_',num2str(evtid)],'half duration',':',true); 
    tshift  = tstart + 2*thd;
else
    % When running SPECFEM alone, seismograms begin at -2x the source
    % duration.
    tstart = sf_read_input([prjDir,'/source_files/CMTSOLUTION_',num2str(evtid)],'half duration',':',true);
    tstart = -2*tstart;
    tshift = 0;
end

% Sampling period/frequency of seismograms
Ts     = sf_read_input([prjDir,'/DATA/Par_file'],'DT','=',true);
fs     = 1/Ts;

% Assign default data window if necessary
if isempty(dataWin)
    nstep   = sf_read_input([prjDir,'/DATA/Par_file'],'NSTEP','=',true);
    dataWin = tstart + [0,nstep*Ts];
end

% Event index
ievt = event.id == evtid;

% Get seismogram file list
theFiles = [];
for ii = 1:length(channel)
    theChannels = dir([prjDir,'/SEIS/CMTSOLUTION_',num2str(evtid),'/*',ext,channel{ii},'_SU']);
    theFiles    = cat(1,theFiles,theChannels);
end

if isempty(theFiles)
    warning('No data found matching criteria');
end

% TO ADD: Restrict number of traces that can be loaded at once

% To ADD: User input sampling frequency

% Check sampling frequency
if abs(1 - (fs/round(fs))) > 0.0001
    warning(['Non integer sampling frequencies not yet supported. Re-sampling from ',num2str(fs),' Hz to ',num2str(round(fs)),' Hz.']);
    fs = round(fs);
    
    % TO ADD: Resampling for non-integer sampling frequencies
else
    % Just forces fs to be an integer
    fs = round(fs);
end

% Construct data structure
% Event parameters
data.event  = evtid;
% Common trace parameters
data.fs     = fs; % Trace sampling frequency
data.ntrace = length(station.id);
data.nchan  = length(channel);
data.nsamp  = 1 + diff(dataWin)*data.fs;
data.t      = dataWin(1) + linspace(0,data.nsamp-1,data.nsamp)'./data.fs;
% Receiver parameters
data.station = station.id;
data.channel = channel;
data.seis    = zeros(data.ntrace,data.nchan,data.nsamp);
% Station range and back-azimuth
if strcmp(tf_coupled,'.true.')
    [data.delta,data.baz] = distance(station.latitude,station.longitude,event.latitude(ievt),event.longitude(ievt));
else
    % If not a coupled simulation, then assume sources are inside the model
    % and source lat/lon corresponds to y/x
    data.delta = sqrt((station.x-(event.longitude(ievt)/1000)).^2 + (station.y-(event.latitude(ievt)/1000)).^2);
    data.baz   = atand((station.x-(event.longitude(ievt)/1000))./(station.y-(event.latitude(ievt)/1000)));
end

% Interpolate travel-times from table
data.phase     = '';
data.tt1D      = zeros(data.ntrace,1);
data.rayP      = zeros(data.ntrace,1);
data.incidence = zeros(data.ntrace,1);
if tf_moveout
    data.phase = aPhase;
    for ii = 1:length(data.delta)
        [tt,rayP,inc] = taup_time(aModel,aPhase,data.delta(ii),event.depth(ievt),0);
        data.tt1D(ii)      = tt;
        data.rayP(ii)      = rayP;
        data.incidence(ii) = inc;
    end
end

% Load and process seismic data
nmissed = 0; % Tracks number of missed stations
stlist  = {};
for ii = 1:length(theFiles)
    % Load the data
    [seis,H] = ReadSu([theFiles(ii).folder,'/',theFiles(ii).name]);
    
    % Identify the channel we loaded (column in j-index in seis array)
    chn = strsplit(theFiles(ii).name,'_');
    chn = chn{2}(end);
    col = strcmp(data.channel,chn);
    
    % Check that we identified a channel
    if sum(col) == 0
        error(['Bad channel. Channel ',chn,' does not match input channel names.']);
    elseif sum(col) > 1
        error(['Bad channel. Channel ',chn,' is not unique.']);
    end
    
    % Loop over traces
    for jj = 1:length(H)
        % Processing j-th trace
        sj = seis(:,jj);
        fs = 1/(H(jj).dt/1000); % SU stores sample interval in milliseconds?
        tj = tstart + linspace(0,H(jj).ns-1,H(jj).ns)'./fs; % Time vector for jth trace
        
        % Identify station based on location
        row = (station.x == H(jj).GroupX/1000) & (station.y == H(jj).GroupY/1000);
        stlist = cat(1,stlist,station.id(row));
        if sum(row) > 1
            error('Multiple stations match identification criteria');
        elseif sum(row) == 0
            nmissed = nmissed + 1;
        else
            
            % Check for consistent sample frequencies
            if abs(1 - (fs/data.fs)) > eps
                % Apply taper to trace ends. This prevents edge artifact when resampling.
                w  = get_window(H(jj).ns,fs,ttap);
                sj = sj.*w(:);
                
                % Resamples with polyphase antialiassing filter
                sj = resample(sj,tj,data.fs);
            end
            
            % Interpolate traces to common time vector
            % sj = interp1(tshift+tj-data.tt1D(row),sj,data.t,'linear',0);
            sj = interp1(tj-data.tt1D(row),sj,data.t,'linear',0);
            
            % Apply a final taper on the trace
            w  = get_window(data.nsamp,data.fs,ttap);
            sj = sj.*w(:);
            
            % Apply filter?
            if ~isempty(filtopts)
                sj = ButterFilter(sj,filtopts(1:2),filtopts(3),data.fs,logical(filtopts(4)));
            end
            
            % Store result
            data.seis(row,col,:) = sj(:);
        end
    end
end
stlist = unique(stlist);

% Rotate seismograms to TRZ?
if tf_rotate && (data.nchan > 1)
    for ii = 1:data.ntrace
        F = [squeeze(data.seis(ii,1,:))'; squeeze(data.seis(ii,2,:))'; zeros(1,data.nsamp)];
        % Rotate into transverse and radial components assuming channel
        % index 1 is east and channel index 2 is north. Note that baz is
        % measured positive clockwise from north while rotz function
        % assumes CCW positive rotations.
        F = rotz(data.baz(ii))*F;
        data.seis(ii,1,:) = F(1,:)'; % Transverse
        data.seis(ii,2,:) = F(2,:)'; % Radial
    end
end

% Did we miss any stations?
% For some reason, every SU file contains an extra trace
nmissed = nmissed - length(theFiles);
if nmissed > 0
    warning(['Seismograms for ',num2str(nmissed),' stations was not found.']);
end
