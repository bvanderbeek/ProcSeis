function data = loadSeis_commonEvent(dataDir,evtid,event,station,sampFreq,...
    timeWindow,refModel,aPhase,corners,order,tf_zerophase)

% Assumes a specific directory structure and file names
% 1) All seismograms are in directories named 'dataDir/EVT#' where # is the
%    event identifier
% 2) All seismograms are miniSEED files named network_station_channel.mseed

% EXPLICITLY ASSUME 3-COMPONENT RECORDS
nchan = 3;

% Parse filenames into network, station, and channel names
theFiles = dir([dataDir,'/EVT',num2str(evtid),'/*.mseed']);
if isempty(theFiles)
    error(['No data found in: ',dataDir,'/EVT',num2str(evtid),'/*.mseed']);
else
    NW_ST_CH = cell(length(theFiles),3);
    for ii = 1:length(theFiles)
        ID = strsplit(theFiles(ii).name,'_');
        NW_ST_CH{ii,1} = ID{1};
        NW_ST_CH{ii,2} = ID{2};
        ID = strsplit(ID{3},'.mseed');
        NW_ST_CH{ii,3} = ID{1};
    end
end

% Initialize seismogram data structure.
% Event and stations for this data structure
data.event.id         = evtid;
data.event.originTime = datenum(event.originDate(event.id == evtid)); % DATENUM
data.event.latitude   = event.latitude(event.id == evtid);
data.event.longitude  = event.longitude(event.id == evtid);
data.event.depth      = event.depth(event.id == evtid);
% Station list
data.network = {};
data.station = unique(NW_ST_CH(:,2));
% Common time vector for all traces
data.nsta  = length(data.station);
data.nsamp = 1 + round(sampFreq*timeWindow);
data.nchan = nchan;
data.fs    = sampFreq;

% Compute range and back azimuth of event
[~,ista] = ismember(data.station,station.name);
[data.delta,data.baz] = distance(station.latitude(ista),station.longitude(ista),...
    data.event.latitude,data.event.longitude);

% Predicted travel-times
TT = zeros(data.nsta,1);
for ii = 1:data.nsta
    % TauP predicted travel-time
    TT(ii) = taup_time(refModel,aPhase,data.delta(ii),data.event.depth,0);
end

% Check for bad arrivals
ibad = isnan(TT);
if sum(ibad) > 0
    warning(['No TauP prediction returned for ',num2str(sum(ibad)),' arrivals; interpolating.']);
    TT(ibad) = interp1(data.delta(~ibad),TT(~ibad),data.delta(ibad),'linear','extrap');
end
% Final check
if any(isnan(TT)) || any(isinf(TT))
    error('Bad predictions!');
end

% Define absolute reference time (DATENUM) for all traces.
data.refTime = data.event.originTime + (TT./(60*60*24));
% Define relative time vector (SECONDS) common to all traces. Note that
% data.t == 0 corresponds to data.refTime.
data.t = (1/data.fs)*linspace(0,data.nsamp-1,data.nsamp) - (timeWindow/2);
% Quick check sampling frequency
if data.fs*max(abs(diff(data.t) - (1/data.fs))) > 1e-6
    error('Time vector is inconsistent with sampling frequency!');
end

% Allocate seismogram array
data.seis = zeros(data.nsta,data.nsamp,data.nchan);
% Loop to read in seismograms
for ii = 1:data.nsta
    % Pointer to file name
    pfl = strcmp(NW_ST_CH(:,2),data.station{ii});
    % Pointer to station
    pst = strcmp(station.name,data.station{ii});
    
    % Identify seismogram files
    ntw  = NW_ST_CH(pfl,1);
    stn  = NW_ST_CH(pfl,2);
    chn  = NW_ST_CH(pfl,3);
    
    % Store network
    if length(unique(ntw)) == 1
        data.network = cat(1,data.network,ntw(1));
    else
        error('Non-unique network name!');
    end
    
    % Initialize channel orientation for current station
    azm = [];
    % Loop over channels present
    for ich = 1:length(chn)
        % Identify channel index
        [~,n] = ismember(chn{ich},station.channel);
        if (n == 1) || (n == 4) || (n == 5) || (n == 8)
            % Component 1 or N. Most network conventions use channel 1 or N
            % as the first component. Here we are reversing this order to
            % be consistent with math.
            kk  = 2;
            % Only this channel defines the instrument orientation
            azm = station.channel_azimuth(pst,n);
        elseif (n == 2) || (n == 3) || (n == 6) || (n == 7)
            % Channel 2 or E
            kk = 1;
        else
            % Channel Z
            if any(strcmp(chn{ich},{'BHZ','HHZ'}))
                kk = 3;
            else
                error(['Unknown channel: ',chn{ich}]);
            end
        end
        
        % Read miniSEED data
        [S,~] = rdmseed([dataDir,'/EVT',num2str(evtid),'/',...
            ntw{ich},'_',stn{ich},'_',chn{ich},'.mseed']);
        
        % Concatonate data blocks
        t  = double(cat(1,S.t)); t = t(:)'; % DATENUM
        d  = double(cat(1,S.d)); d = d(:)';
        fs = double(cat(1,S.SampleRate)); fs = fs(:)';
        if (range(fs)/mean(fs)) > (1e-10)
            error('Non-unique sampling frequencies!');
        else
            fs = mean(fs);
        end
        
        % Time with respect to predicted arrival converted to seconds
        t = (t - data.refTime(ii))*24*60*60;
        
        % Resample if necessary
        if ~(fs == data.fs)
            d = resample(d,data.fs,fs);
            t = t(1) + (1/data.fs)*linspace(0,length(d)-1,length(d));
        end
        
        % Bandpass filter
        d = ButterFilter(d,corners,order,data.fs,tf_zerophase);
        
        % Interpolate to common time vector
        d = interp1(t,d,data.t,'linear',0);
        
        % Store seismograms
        data.seis(ii,:,kk) = d;
    end
    
    if isempty(azm) || isnan(azm)
        warning('Bad orientation!');
        keyboard;
    end
    
    % Rotate Seismograms
    
    S12 = squeeze(data.seis(ii,:,1:2))';
    if ~(azm == 0)
        % Rotate from 1,2 to E,N. This is a clockwise rotation and thus the
        % sign of azm is flipped. This is because we are rotating from 1,2-
        % to E,N-coordinates but azm is the measured in geographic
        % coordinates clockwise of north (draw it and it makes sense).
        R   = [cosd(-azm), -sind(-azm); sind(-azm), cosd(-azm)]; % CCW negative convention
        S12 = R*S12;
    end
    % Rotate from E,N- to TRZ-coordinates. This is a counterclockwise
    % rotation by baz. This is because we are rotating from E,N- to
    % TRZ-coordinates but baz is measured in the geographic coordinates
    % clockwise of N.
    R   = [cosd(data.baz(ii)), -sind(data.baz(ii));...
           sind(data.baz(ii)),  cosd(data.baz(ii))]; % CCW negative convention
    S12 = R*S12;
    % Update seismograms
    data.seis(ii,:,1) = S12(1,:)';
    data.seis(ii,:,2) = S12(2,:)';
    
    fprintf(['Loaded ',num2str(ii),' of ',num2str(data.nsta),' stations. \n']);
end