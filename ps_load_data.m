function handles = ps_load_data(handles)

DataDir = ps_get_value(handles,'DataDir');
Ext     = ps_get_value(handles,'Ext');

if strcmp(Ext,'mseed')
    % Define event data
    handles.event.id         = zeros(nsrc,1); % Event id (sequentially numbered from file name)
    handles.event.origintime = zeros(nsrc,1); % Event origin time (all SPECFEM sources start at 0 s)
    handles.event.latitude   = zeros(nsrc,1); % Event latitude (dd)
    handles.event.longitude  = zeros(nsrc,1); % Event longitude (dd)
    handles.event.depth      = zeros(nsrc,1); % Event depth measured with respect to model surface (km)
    handles.event.Mw         = zeros(nsrc,1); % Event moment magnitude
    
    % Define station date
    handles.station.id        = cell(nrec,1); % Station ID
    handles.station.longitude = zeros(nrec,1); % (dec. deg)
    handles.station.latitude  = zeros(nrec,1); % (dec. deg)
    handles.station.elevation = zeros(nrec,1); % Receiver elevation (km)...this is the SPECFEM z-coordinate measured with respect to model surface
    handles.station.depth     = zeros(nrec,1); % Receiver depth with respect to model surface (km)
    
    % Define seismic data
    % Extract some information from GUI menues
    [ChannelList,MinTime,MaxTime,aModel,aPhase] = ps_get_value(handles,...
        'ChannelList','MinTime','MaxTime','RefModel','RefPhase');
    handles.SeisDat.event   = 0;
    handles.SeisDat.channel = ChannelList; % Cell array of channel names that are loaded
    handles.SeisDat.fs      = Fs; % Sampling frequency that is common to all traces
    handles.SeisDat.ntrace  = NTRACES; % Number of seismic traces loaded
    handles.SeisDat.nchan   = NCHAN; % Number of channels in each trace
    handles.SeisDat.nsamp   = 1 + round(Fs*(MaxTime - MinTime)); % Number of samples in each trace
    handles.SeisDat.t       = MinTime + (1/Fs)*linspace(0,handles.SeisDat.nsamp - 1,handles.SeisDat.nsamp)'; % Reduced time vector for all traces
    handles.SeisDat.phase   = repmat({aPhase},handles.SeisDat.ntrace); % Seismic phase reference for reduced time
    % Allocate variables to be determined
    handles.SeisDat.station   = cell(handles.SeisDat.ntrace,1); % Station corresponding to each trace
    handles.SeisDat.delta     = zeros(handles.SeisDat.ntrace,1); % Source-receiver range for each trace (deg.)
    handles.SeisDat.baz       = zeros(handles.SeisDat.ntrace,1); % Source-receiver back-azimuth for each trace (deg.)
    handles.SeisDat.tt1D      = zeros(handles.SeisDat.ntrace,1); % TauP predicted time
    handles.SeisDat.rayP      = zeros(handles.SeisDat.ntrace,1); % TauP predicted ray parameter
    handles.SeisDat.incidence = zeros(handles.SeisDat.ntrace,1); % TauP predicted incidence angle
    handles.SeisDat.seis      = zeros(handles.SeisDat.ntrace,handles.SeisDat.nchan,handles.SeisDat.nsamp); % Seismogram array
    
elseif strcmp(Ext,'mat')
    % Load events structure
    theEvents = dir([DataDir,'/Events_*']);
    if length(theEvents) == 1
        load([DataDir,'/',theEvents.name],'Events');
    elseif length(theEvents) > 1
        error('Multiple event files found.');
    else
        error('Cannot locate event file.');
    end
    % Pre-allocate event structure
    nsrc = length(Events);
    handles.event.id         = zeros(nsrc,1); % Event id (sequentially numbered from file name)
    handles.event.origintime = zeros(nsrc,1); % Event origin time (all SPECFEM sources start at 0 s)
    handles.event.latitude   = zeros(nsrc,1); % Event latitude (dd)
    handles.event.longitude  = zeros(nsrc,1); % Event longitude (dd)
    handles.event.depth      = zeros(nsrc,1); % Event depth measured with respect to model surface (km)
    handles.event.Mw         = zeros(nsrc,1); % Event moment magnitude
    for ii = 1:length(Events)
        % Extract ID
        EID = strsplit(Events(ii).PublicId,'=');
        EID = str2double(EID{end});
        % Fill structure
        handles.event.id(ii)         = EID;
        handles.event.origintime(ii) = 60*60*24*datenum(Events(ii).PreferredTime,'yyyy-mm-dd HH:MM:SS.FFF');
        handles.event.latitude(ii)   = Events(ii).PreferredLatitude;
        handles.event.longitude(ii)  = Events(ii).PreferredLongitude;
        handles.event.depth(ii)      = Events(ii).PreferredDepth;
        handles.event.Mw(ii)         = Events(ii).PreferredMagnitude.Value;
    end
    
    % Load stations structure
    theStations = dir([DataDir,'/Stations_*']);
    if length(theStations) == 1
        load([DataDir,'/',theStations.name],'Stations');
    elseif length(theStations) > 1
        error('Multiple station files found.');
    else
        error('Cannot locate station file.');
    end
    % Pre-allocate station structure
    nrec = length(Stations);
    handles.station.id        = cell(nrec,1); % Station ID
    handles.station.longitude = zeros(nrec,1); % (dec. deg)
    handles.station.latitude  = zeros(nrec,1); % (dec. deg)
    handles.station.elevation = zeros(nrec,1); % Receiver elevation (km)...this is the SPECFEM z-coordinate measured with respect to model surface
    handles.station.depth     = zeros(nrec,1); % Receiver depth with respect to model surface (km)
    % Fill stations
    for ii = 1:nrec
        handles.station.id{ii}        = [Stations(ii).NetworkCode,'_',Stations(ii).StationCode];
        handles.station.longitude(ii) = Stations(ii).Longitude;
        handles.station.latitude(ii)  = Stations(ii).Latitude;
        handles.station.elevation(ii) = Stations(ii).Elevation/1000;
    end
    
    % Load the seismograms
    theTraces = dir([DataDir,'/IFTrace_*']);
    if isempty(theTraces)
        error('No traces found.');
    end
    load([DataDir,'/',theTraces.name],'Traces'); % FIX ME!
    
    % Determine event from filename
    eind = strsplit(theTraces.name,'_');
    eind = strsplit(eind{end},'.mat');
    eind = str2double(eind{1});
    % Determine highest sampling frequency
    Fs = [Traces{:}];
    Fs = max([Fs(:).sampleRate]);
    % Extract some information
    [ChannelList,MinTime,MaxTime,aModel,aPhase] = ps_get_value(handles,...
        'ChannelList','MinTime','MaxTime','RefModel','RefPhase');
    
    % Parameters common to all traces
    handles.SeisDat.event   = handles.event.id(eind);
    handles.SeisDat.channel = ChannelList;
    handles.SeisDat.fs      = Fs; 
    handles.SeisDat.ntrace  = length(Traces);
    handles.SeisDat.nchan   = length(ChannelList);
    handles.SeisDat.nsamp   = 1 + round(Fs*(MaxTime - MinTime));
    handles.SeisDat.t       = MinTime + (1/Fs)*linspace(0,handles.SeisDat.nsamp - 1,handles.SeisDat.nsamp)';
    handles.SeisDat.phase   = repmat({aPhase},handles.SeisDat.ntrace);
    % Allocate variables to be determined
    handles.SeisDat.station   = cell(handles.SeisDat.ntrace,1);
    handles.SeisDat.delta     = zeros(handles.SeisDat.ntrace,1);
    handles.SeisDat.baz       = zeros(handles.SeisDat.ntrace,1);
    handles.SeisDat.tt1D      = zeros(handles.SeisDat.ntrace,1);
    handles.SeisDat.rayP      = zeros(handles.SeisDat.ntrace,1);
    handles.SeisDat.incidence = zeros(handles.SeisDat.ntrace,1);
    handles.SeisDat.seis      = zeros(handles.SeisDat.ntrace,handles.SeisDat.nchan,handles.SeisDat.nsamp);
    % Fill fields
    for ii = 1:handles.SeisDat.ntrace
        aTrace = Traces{ii};
        for jj = 1:length(aTrace)
            nsj = aTrace(jj).sampleCount;
            fsj = aTrace(jj).sampleRate;
            tij = 60*60*24*aTrace(jj).startTime;
            
            % Define absolute time vector
            tt  = tij + (1/fsj)*linspace(0,nsj-1,nsj)';
            if (abs(tt(1) - tij) > (1e-4)) ||...
                    (abs(tt(end) - 60*60*24*aTrace(jj).endTime) > (1e-4))
                error('Bad time vector.');
            end
            % Remove origin time
            tt = tt - handles.event.origintime(eind);
            
            % Distance and back-azimuth
            [D,BAZ] = distance(aTrace(jj).latitude,aTrace(jj).longitude,...
                handles.event.latitude(eind),handles.event.longitude(eind));
            handles.SeisDat.delta(ii) = D;
            handles.SeisDat.baz(ii)   = BAZ;
            
            % Compute travel-time
            [tt1D,rayP,inc] = taup_time(aModel,aPhase,D,handles.event.depth(eind),0);
            handles.SeisDat.tt1D(ii)      = tt1D;
            handles.SeisDat.rayP(ii)      = rayP;
            handles.SeisDat.incidence(ii) = inc;
            
            % Interpolate traces to common time vector
            s = interp1(tt - tt1D,aTrace(jj).data,handles.SeisDat.t,'linear',0);
            s = s./aTrace(jj).sensitivity;
            
            % Filter
            s = ButterFilter(s,[1/100,1/10],2,Fs,true);
            
            % Store data
            handles.SeisDat.seis(ii,jj,:) = s;
        end
    end
else
    handles.station = sf_read_stations(DataDir);
    handles.event   = sf_read_events(DataDir);
    handles         = ps_load_SFseis(handles);
end