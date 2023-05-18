function data = sf_load_ascii_seis(station,event,evtid,channel,prjDir,ext,dataWin,varargin)


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
    error('Need to input aPhase and a 1D model name to use moveout correction.');
end

% Get seismogram start times. This is defined differently for an AxiSEM
% coupled simulation.
tf_coupled = sf_read_input([prjDir,'/DATA/Par_file'],'COUPLE_WITH_INJECTION_TECHNIQUE','=',false);
ctype      = sf_read_input([prjDir,'/DATA/Par_file'],'INJECTION_TECHNIQUE_TYPE','=',true);
if strcmp(tf_coupled,'.true.') && (ctype == 2)
    tstart  = dlmread([prjDir,'/','/source_files/reformat_CMTSOLUTION_',num2str(evtid),'.par']);
    tstart  = tstart(2,1);
    % SPECFEM seismograms start at -2x the source duration
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

% Check sampling frequency
if abs(1 - (fs/round(fs))) > 0.0001
    warning(['Non integer sampling frequencies not yet supported. Re-sampling from ',num2str(fs),' Hz to ',num2str(round(fs)),' Hz.']);
    fs = round(fs);
    
    % TO ADD: Resampling for non-integer sampling frequencies
else
    % Just forces fs to be an integer
    fs = round(fs);
end

% Event index
ievt = event.id == evtid;

% Get seismogram file list
theFiles = [];
ntrace   = 0;
station_list = {};
for ii = 1:length(channel)
    theChannels = dir([prjDir,'/SEIS/CMTSOLUTION_',num2str(evtid),'/*',channel{ii},'.',ext]);
    % Remove channels not in station array
    irmv = false(length(theChannels),1);
    stid = cell(length(theChannels),1);
    for jj = 1:length(theChannels)
        stinfo   = strsplit(theChannels(jj).name,'.');
        stid{jj} = [stinfo{1},stinfo{2}];
        if sum(strcmp(station.id,stid{jj})) == 0
            irmv(jj) = true;
        end
    end
    station_list      = cat(1,station_list,stid(~irmv));
    theChannels(irmv) = [];
    ntrace            = max(ntrace,length(theChannels));
    theFiles          = cat(1,theFiles,theChannels);
end
% Remove stations without seismograms
irmv = ~ismember(station.id,station_list);
station.id(irmv)        = [];
station.x(irmv)         = [];
station.y(irmv)         = [];
station.elevation(irmv) = [];
station.depth(irmv)     = [];
station.longitude(irmv) = [];
station.latitude(irmv)  = [];

if isempty(theFiles)
    warning('No data found matching criteria');
end

% TO ADD: Restrict number of traces that can be loaded at once

% To ADD: User input sampling frequency

% Construct data structure
% Event parameters
data.event  = evtid;
% Common trace parameters
data.fs     = fs; % Trace sampling frequency
data.ntrace = ntrace;
data.nchan  = length(channel);
data.nsamp  = diff(dataWin)*data.fs + 1;
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

% Use TauP to calculate 1D travel-times
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
for ii = 1:length(theFiles)
    % Derive station and channel from file name
    stinfo = strsplit(theFiles(ii).name,'.');
    stid   = [stinfo{1},stinfo{2}];
    chid   = stinfo{3};
    
    % Index of trace in data array
    row = strcmp(data.station,stid);
    col = strcmp(data.channel,chid);
    
    % Load the data
    seis = dlmread([theFiles(ii).folder,'/',theFiles(ii).name]);
    
    % Check for consistent sample frequencies
    fs = 1/mean(unique(diff(seis(:,1))));
    if abs(1 - (fs/data.fs)) > 0.0001
        % Apply taper to trace ends. This prevents edge artifact when resampling.
        w = get_window(length(seis(:,1)),fs,ttap);
        seis(:,2) = seis(:,2).*w(:);
        
        % Resamples with polyphase antialiassing filter
        seis(:,2) = resample(seis(:,2),seis(:,1),data.fs);
    end
    
    % Interpolate traces to common time vector
    seis = interp1(tshift+seis(:,1)-data.tt1D(row),seis(:,2),data.t,'linear',0);
    
    % Apply a final taper on the trace
    w    = get_window(data.nsamp,data.fs,ttap);
    seis = seis.*w(:);
    
    % Apply filter?
    if ~isempty(filtopts)
        seis = ButterFilter(seis,filtopts(1:2),filtopts(3),data.fs,logical(filtopts(4)));
    end
    
    % Store result
    data.seis(row,col,:) = seis(:);
end

% Rotate seismograms to TRZ?
if isstruct(tf_rotate) %#ok Is defined by eval statement at beginning
    R                = tf_rotate;
    tf_rotate        = true;
    [ll,R.ista]      = ismember(R.station,data.station);
    data.orientation = zeros(data.ntrace,1);
    data.orientation(R.ista) = R.orientation(ll);
elseif tf_rotate
    R.station     = data.station;
    R.orientation = data.baz;
    R.ista        = (1:data.ntrace)';
end
if tf_rotate && (data.nchan > 1)
    for ii = 1:data.ntrace
        F = [squeeze(data.seis(ii,1,:))'; squeeze(data.seis(ii,2,:))'; zeros(1,data.nsamp)];
        % Rotate by +baz since baz is measured CW from N (rotz angle is CCW
        % from +x). Channel 1 is the transverse and channel 2 is radial.
        F = rotz(R.orientation(R.ista(ii)))*F;
        data.seis(ii,1,:) = F(1,:)';
        data.seis(ii,2,:) = F(2,:)';
    end
end

% OLD rotation option
% % Rotate seismograms to TRZ?
% if tf_rotate && (data.nchan > 1)
%     for ii = 1:data.ntrace
%         F = [squeeze(data.seis(ii,1,:))'; squeeze(data.seis(ii,2,:))'; zeros(1,data.nsamp)];
%         % Rotate by +baz since baz is measured CW from N (rotz angle is CCW
%         % from +x). Channel 1 is the transverse and channel 2 is radial.
%         F = rotz(data.baz(ii))*F;
%         data.seis(ii,1,:) = F(1,:)';
%         data.seis(ii,2,:) = F(2,:)';
%     end
% end
