function handles = ps_load_mseed_data(handles)
% Hard-coded sampling frequency!!!
sampFreq = 40;

% Parse data loading parameters
data_dir = ps_get_value(handles,'DataDir'); % Directory with data
ievt = str2double(handles.DTxtEind.String); % Event index (NOT ID) to load
dataWin = [str2double(handles.DTxtMinTime.String),str2double(handles.DTxtMaxTime.String)]; % Time window
fmin    = eval(handles.DTxtMinFreq.String); % Min filter frequency
fmax    = eval(handles.DTxtMaxFreq.String); % Max filter frequency
order   = str2double(handles.DTxtOrder.String); % Filter order
tf_zero = handles.ChkZeroPhase.Value; % Zero-phase filter?

% Load pre-defined event and station structures
load([data_dir,'/events.mat'],'events');
load([data_dir,'/stations.mat'],'stations');

% Store event and station data
handles.event = events;
handles.station = stations;
% Re-define station 'names' to station 'ids'
handles.station.id = stations.name;

% Extract some information from GUI menues
[refModel,aPhase] = ps_get_value(handles,'RefModel','RefPhase');

handles.SeisDat = loadSeis_commonEvent(data_dir,handles.event.id(ievt),handles.event,handles.station,sampFreq,...
    diff(dataWin),refModel,aPhase,[fmin,fmax],order,tf_zero);


handles.SeisDat.event   = handles.event.id(ievt);
handles.SeisDat.channel = {'E', 'N', 'Z'}; % Cell array of channel names that are loaded
handles.SeisDat.fs      = sampFreq; % Sampling frequency that is common to all traces
handles.SeisDat.ntrace  = handles.SeisDat.nsta; % Number of seismic traces loaded
handles.SeisDat.phase   = repmat({aPhase},handles.SeisDat.ntrace); % Seismic phase reference for reduced time
% Allocate variables to be determined
handles.SeisDat.tt1D    = (handles.SeisDat.refTime - datenum(handles.event.originDate(ievt)))*60*60*24;
handles.SeisDat.seis    = permute(handles.SeisDat.seis, [1,3,2]);
handles.SeisDat.t = handles.SeisDat.t(:);

handles.MenuCoord.String = 'ENZ';
