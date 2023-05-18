function event = sf_read_events(prjDir)

% LIMITATIONS:
% + Assumes event locations are in geographic coordinates. If this is not
%   the case, then longitude is x- and latitude is y-coordinate
%

% Index event files
theFiles = dir([prjDir,'/source_files/CMT*']);
nsrc     = length(theFiles);

% Pre-allocate event structure
event.id         = zeros(nsrc,1); % Event id (sequentially numbered from file name)
event.latitude   = zeros(nsrc,1); % Event latitude (dd)
event.longitude  = zeros(nsrc,1); % Event longitude (dd)
event.depth      = zeros(nsrc,1); % Event depth measured with respect to model surface (km)
event.origintime = zeros(nsrc,1); % Event origin time (all SPECFEM sources start at 0 s)
event.Mo         = zeros(nsrc,1);
event.Mw         = zeros(nsrc,1); % Event moment magnitude
for ii = 1:length(theFiles)
    id  = strsplit(theFiles(ii).name,'_');
    id  = str2double(id{2});
    event.id(ii)        = id;
    event.latitude(ii)  = sf_read_input([theFiles(ii).folder,'/',theFiles(ii).name],'latorUTM',':',true);
    event.longitude(ii) = sf_read_input([theFiles(ii).folder,'/',theFiles(ii).name],'longorUTM',':',true);
    event.depth(ii)     = sf_read_input([theFiles(ii).folder,'/',theFiles(ii).name],'depth',':',true);
    
    % Magnitude calculation from tensor components (see SPECFEM manual on
    % CMTSOLUTION file)
    Mrr = sf_read_input([theFiles(ii).folder,'/',theFiles(ii).name],'Mrr',':',true);
    Mtt = sf_read_input([theFiles(ii).folder,'/',theFiles(ii).name],'Mtt',':',true);
    Mpp = sf_read_input([theFiles(ii).folder,'/',theFiles(ii).name],'Mpp',':',true);
    Mrt = sf_read_input([theFiles(ii).folder,'/',theFiles(ii).name],'Mrt',':',true);
    Mrp = sf_read_input([theFiles(ii).folder,'/',theFiles(ii).name],'Mrp',':',true);
    Mtp = sf_read_input([theFiles(ii).folder,'/',theFiles(ii).name],'Mtp',':',true);
    Mo  = sqrt(Mrr^2 + Mtt^2 + Mpp^2 + 2*Mrt^2 + 2*Mrp^2 + 2*Mtp^2)/sqrt(2);
    event.Mo(ii) = Mo/(10^7); % Moment in N*m
    event.Mw(ii) = 2*(log10(Mo) - 16.1)/3;
    
%     % Event magnitude from header (10 index is body wave magnitude; 11
%     % index is surface wave magnitude)
%     H   = strsplit(sf_read_input([theFiles(ii).folder,'/',theFiles(ii).name],'PDE','PDE',false),' '); % Event header
%     event.magnitude(ii) = str2double(H{10});
end
