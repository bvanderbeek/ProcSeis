function station = sf_read_stations(prjDir,varargin)

% LIMITATIONS:
% + Assumes stations locations are in cartesian coordinates
% + Does not account for a rotated specfem model
%


% Geographic center of SPECFEM model
if isempty(varargin)
    if isfile([prjDir,'/MESH/ParFileMeshChunk'])
        params       = sf_read_lines([prjDir,'/MESH/ParFileMeshChunk']);
        lon0lat0azim = str2num(params{2}); %#ok <is an array>
    else
        warning('Could not find /MESH/ParFileMeshChunk. Assuming origin of model at 0 Lat, 0 Lon.');
        lon0lat0azim = [0 0 0];
    end
elseif length(varargin) == 1
    lon0lat0azim = varargin{1};
else
    error('Incorrect number of variable input arguments.');
end

if ~(lon0lat0azim(3) == 0)
    error('Rotation of model domain is not yet supported');
end

% Receiver parameters
params = sf_read_lines([prjDir,'/DATA/STATIONS']);
nrec   = length(params);
% Pre-allocate station structure
station.id        = cell(nrec,1); % Station ID
station.x         = zeros(nrec,1); % Cartesian x-coordinate (km)
station.y         = zeros(nrec,1); % Cartesian y-coordinate (km)
station.elevation = zeros(nrec,1); % Receiver elevation (km)...this is the SPECFEM z-coordinate measured with respect to model surface
station.depth     = zeros(nrec,1); % Receiver depth with respect to model surface (km)
for ii = 1:nrec
    % Interpret station info
    stinfo   = strsplit(strtrim(params{ii}));
    stid     = stinfo{1}; % Station name
    stnw     = stinfo{2}; % Network code
    station.id{ii}        = [stnw,stid];
    station.x(ii)         = str2double(stinfo{4})/1000;
    station.y(ii)         = str2double(stinfo{3})/1000;
    station.elevation(ii) = str2double(stinfo{5})/1000;
    station.depth(ii)     = str2double(stinfo{6})/1000;
end
% Derive geographic coordinates
[station.longitude,station.latitude] = cart2geo(station.x,station.y,station.elevation,lon0lat0azim(1),lon0lat0azim(2),true);
