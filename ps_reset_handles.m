function handles = ps_reset_handles(handles)

% Define default parameters
if ~isfield(handles,'tf_reset_defaults')
    handles.tf_reset_defaults = true;
end

if handles.tf_reset_defaults
    % Reset all analysis parameters (e.g., when the program is first
    % started)
    
    % Define default parameters
    run('set_defaults_ProcSeis.m');
    
    % Add required matlab tools to path
    procseis_path = mfilename('fullpath');
    procseis_path = strsplit(procseis_path,'ps_reset_handles');
    addpath(genpath([procseis_path{1},'contrib']));
    % Define TauP environment variable
    setenv('TAUPJAR',[procseis_path{1},'contrib/TauP-2.4.5.jar']);
    
    % Set GUI Default Values
    Name  = fieldnames(Defaults);
    Value = struct2cell(Defaults);
    % Reset values
    handles = ps_set_value(handles,Name,Value);
    
    % Set some additional values
    handles = ps_set_value(handles,'EventIndex',1);
    
    % Load Data
    handles = ps_load_data(handles);
    
    % Display Panel
    % Dynamic Text
    handles = ps_set_value(handles,{'XMin','XMax'},{Defaults.MinTime,Defaults.MaxTime});
    
    % Track that defaults were already initialized
    handles.tf_reset_defaults = false;
else
    % Do not reset all parameters (e.g. when loading new seismograms)
    
    % Load Data
    handles = ps_load_data(handles);
    
    % Reset trace indexing
    NTraces = ps_get_value(handles,'TraceIndices');
    NTraces = min(length(NTraces),handles.SeisDat.ntrace);
    handles = ps_set_value(handles,'TraceIndices',1:NTraces);
end

% Display total number of events
handles.STxtEInd.String = ['Event Index (of ',num2str(length(handles.event.id)),')'];

% Pop-up Menus
% Trace sorting options
handles = ps_set_value(handles,'SortName',{'traceid','delta','baz'}');
handles = ps_set_value(handles,'SortIndex',1);

% Orientation options
% Custom orientations
if isfield(handles.station,'rotation')
    xyz  = fieldnames(handles.station.rotation);
    ienz = strcmp(xyz,'ENZ');
    xyz  = xyz(~ienz);
else
    xyz = {};
end
% Define available coordinate systems based on number of channels present
if handles.SeisDat.nchan == 1
    CoordName = {'ENZ'};
elseif handles.SeisDat.nchan == 2
    CoordName = cat(1,{'ENZ';'TRZ';'PAZ'},xyz);
elseif handles.SeisDat.nchan == 3
    CoordName = cat(1,{'ENZ';'TRZ';'PAZ';'TQL';'PPD'},xyz);
else
    warning('No rotation options because more than 3 channels loaded.');
    CoordName = {''};
end
handles = ps_set_value(handles,'CoordName',CoordName);
handles = ps_set_value(handles,'CoordIndex',1);

% Track seismometer rotation matrices
if (handles.SeisDat.nchan > 1) && (handles.SeisDat.nchan <= 3)
    % Initialize identity matrices
    Rmat = eye(handles.SeisDat.nchan);
    handles.SeisDat.Rmat = repmat(Rmat(:)',handles.SeisDat.ntrace,1);
end

% Initialize pick sub-structure
handles.picks.t_stack = 0;
handles.picks.DT      = zeros(handles.SeisDat.ntrace,1);
handles.picks.ddt     = zeros(handles.SeisDat.ntrace,1);
handles.picks.CCF     = zeros(handles.SeisDat.ntrace,1);
handles.picks.tf_bad  = false(handles.SeisDat.ntrace,1);
handles.picks.tf_mcc  = false;
handles.picks.SI      = [];
