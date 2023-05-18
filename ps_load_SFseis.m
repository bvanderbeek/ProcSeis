function handles = ps_load_SFseis(handles)
fprintf('\n Loading seismic data...\n');

% Determine seismogram format
tf_bin = sf_read_input([handles.DTxtDataDir.String,'/DATA/Par_File'],...
    'USE_BINARY_FOR_SEISMOGRAMS','=',false);
if strcmp(tf_bin,'.true.')
    handles.tf_bin_seis = true;
elseif strcmp(tf_bin,'.false.')
    handles.tf_bin_seis = false;
else
    error('Problem determining seismogram format.');
end

% Get Event Index
ievt = str2double(handles.DTxtEind.String);

% Data window to load
dataWin = [str2double(handles.DTxtMinTime.String),str2double(handles.DTxtMaxTime.String)];
% Filter parameters
filtopts = [];
if ~isempty(handles.DTxtMinFreq.String) && ~isempty(handles.DTxtMaxFreq.String)...
        && ~isempty(handles.DTxtOrder.String)
    fmin    = eval(handles.DTxtMinFreq.String);
    fmax    = eval(handles.DTxtMaxFreq.String);
    order   = str2double(handles.DTxtOrder.String);
    tf_zero = handles.ChkZeroPhase.Value;
    % Check that valid filter options were returned above
    if ~isempty(fmin) && ~isempty(fmax) && ~isempty(order)
        if (fmin > 0) && (fmax > 0) && ~isnan(fmin) && ~isnan(fmax)...
                && (order > 0) && ~isnan(order)
            filtopts = [fmin,fmax,order,tf_zero];
        end
    end
end
% Channels
chans = strsplit(handles.DTxtChan.String,',');

% Load Seismograms
if handles.tf_bin_seis
    handles.SeisDat = sf_load_bin_seis(handles.station,handles.event,handles.event.id(ievt),...
        chans,handles.DTxtDataDir.String,handles.DTxtExt.String,dataWin,...
        'filtopts',filtopts,'tf_moveout',true,'aPhase',handles.DTxtPhase.String,...
        'aModel',handles.DTxtRefMod.String);
else
    error('Option to load ASCII seismograms not yet implemented.');
end

% Check for additional coordinate systems/rotations
if isfield(handles.station,'rotation')
    chans = fieldnames(handles.station.rotation);
    ienz  = find(strcmp(chans,'ENZ'),1);
    if ~isempty(ienz)
        % Rotate seismograms to East-North-Z using user-define
        % station-specific rotation matrix. The rotation matrix (R) should 
        % be defined such R*[u1;u2;u3] = [ue;un;uz]. In this case, the
        % rotation angle is the counter-clockwise positive angle between x1
        % and true east measured from x1.
        fprintf('\n Applying user-defined rotation to go from XYZ to ENZ coordinates. \n');
        for ii = 1:handles.SeisDat.ntrace
            % Index station
            [~,ista] = ismember(handles.SeisDat.station(ii),handles.station.id);
            % Define rotation matrix
            R = reshape(handles.station.rotation.ENZ(ista,:),3,3);
            R = R(1:handles.SeisDat.nchan,1:handles.SeisDat.nchan);
            % Apply rotation
            handles.SeisDat.seis(ii,:,:) = R*squeeze(handles.SeisDat.seis(ii,:,:));
        end
    end
    % Remove ENZ from additional channel list
    chans(ienz) = [];
else
    chans = {};
end

% Update Menu Options
handles.MenuChan.String = strsplit(handles.DTxtChan.String,',')';
if iscell(handles.MenuChan.String)
    if length(handles.MenuChan.String) == 2
        handles.MenuCoord.String = cat(1,{'ENZ';'TRZ'},chans(:));
    elseif length(handles.MenuChan.String) == 3
        handles.MenuCoord.String = cat(1,{'ENZ';'TRZ';'TQL'},chans(:));
    else
        handles.MenuCoord.String = 'ENZ';
    end
else
    handles.MenuCoord.String = 'ENZ';
end

fprintf(' ...finished loading.\n');