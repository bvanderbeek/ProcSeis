%% WRITE DERIVED FILES
%
% Coupling AxiSEM and SPECFEM3D requires the definition of multiple
% additional input files most of which can be derived from the most basic
% required input files; that is what this script does. Based on existing
% input files this script creates the following:
% 1) MESH/ParFileMeshChunk
% 2) source_files/reformat.par (one for each source)
% 3) source_files/expand_2D_3D.par (one for each source)
% 4) source_files/inparam_source (one for each source, four for each source if using a double-couple)
% 5) A dummy STATION file for AxiSEM (required input but does nothing when
%    coupling with SPECFEM)
%
% There are a couple of baked-in assumptions. These are...
% + I have assumed that the origin of the SPECFEM3D cartesian model is
%   0 latitude and 0 longitude and centered within the grid.
%
%
clear; close all; clc

% Add path to toolbox
addpath(genpath('/Users/bvanderbeek/research/software/sf_toolbox/m-files'));
% Add path to TauP for estimating wavefront arrival
setenv('TAUPJAR','/Users/bvanderbeek/research/software/sf_toolbox/TauP/lib/TauP-2.4.5.jar');
addpath(genpath('/Users/bvanderbeek/research/software/sf_toolbox/TauP/bpv_matlab_toolbox_for_taup'));

% INPUTS
% The main project directory
prjDir  = '/Users/bvanderbeek/research/software/sf_toolbox/examples/Smallbox';
% Cartesian model parameters
R_km   = 6371; % Radius of Earth (km)
azim   = 0; % Azimuth of model domain (deg.); KEEP 0! NOT SURE HOW TO USE OTHERWISE!
lat0   = 42; % Latitude of cartesian model origin (deg.);
lon0   = 15; % Longitude of cartesian model origin (deg.)
% Parameters required for MESH/ParFileMeshChunk
subsur = '.false.'; % If .true., then model domain does not extend to surface; KEEP FALSE!...or make updates

% Model and TauP directory for arrival time prediction
% Begin storing tractions at model edges this many seconds before the tauP 
% predicted arrival time.
a1Dmodel  = 'iasp91'; % 1D Earth model for predicting phase arrivals
aPhase_init = 'P'; % Start storing tractions just before this phase arrives
dt_off      = 25; % Start simulation this many seconds before aPhase_start arrives at bottom edge of model
% Latest phase arrival of interest. Code will check that enough simulation
% time has been defined to capture 'aPhase_fin'
aPhase_fin  = 'S';

%% READ SPECFEM3D PARAMETERS

% The main parameter file
theFile       = [prjDir,'/DATA/Par_file'];
NSTEP         = sf_read_input(theFile,'NSTEP','=',true); % Number of integration time steps
DT            = sf_read_input(theFile,'DT','=',true); % Time step (s)
NPROC         = sf_read_input(theFile,'NPROC','=',true); % Number of processors
LOCAL_PATH    = sf_read_input(theFile,'LOCAL_PATH','=',false); % Path to DATABASES_MPI directory
TRACTION_PATH = sf_read_input(theFile,'TRACTION_PATH','=',false); % Path to where AxiSEM will store tractions


% Mesh parameter files
theFile  = [prjDir,'/DATA/meshfem3D_files/Mesh_Par_file'];
ymin_m   = sf_read_input(theFile,'LATITUDE_MIN','=',true); % in meters!
ymax_m   = sf_read_input(theFile,'LATITUDE_MAX','=',true); % in meters!
xmin_m   = sf_read_input(theFile,'LONGITUDE_MIN','=',true); % in meters!
xmax_m   = sf_read_input(theFile,'LONGITUDE_MAX','=',true); % in meters!
depth_km = sf_read_input(theFile,'DEPTH_BLOCK_KM','=',true); % in *kilometers*!
NX       = sf_read_input(theFile,'NEX_XI','=',true);
NY       = sf_read_input(theFile,'NEX_ETA','=',true);
% Need to read interfaces file for the number of elements in z
theFile = [prjDir,'/DATA/meshfem3D_files/interfaces.dat'];
params  = sf_read_lines(theFile);
Nint    = str2double(params{1}); % First entry is number of interfaces
% Count the number of elements in Z-direction. These are stored in the last
% N-lines of the interfaces.dat file.
NZ = 0;
for ii = (length(params)-Nint+1):length(params)
    NZ = NZ + str2double(params{ii});
end

% Source type expected by AxiSEM
theFile = [prjDir,'/AXISEM_INPUT/inparam_basic'];
SIMULATION_TYPE   = sf_read_input(theFile,'SIMULATION_TYPE','SIMULATION_TYPE',false);

% Number of processors used in AxiSEM calculation
theFile = [prjDir,'/AXISEM_INPUT/inparam_mesh'];
NTHETA  = sf_read_input(theFile,'NTHETA_SLICES','NTHETA_SLICES',true);
NRADIAL = sf_read_input(theFile,'NRADIAL_SLICES','NRADIAL_SLICES',true);

% Load event and station data
event   = sf_read_events(prjDir);
station = sf_read_stations(prjDir,[lon0,lat0,azim]);
nsrc    = length(event.id);

% Check source type
if strcmp(SIMULATION_TYPE,'single')
    SOURCE_TYPE      = {'explosion'};
    HARD_DIR         = {''};
    NSOURCE          = 1;
    SOURCE_AMPLITUDE = event.Mo;
elseif strcmp(SIMULATION_TYPE,'moment')
    SOURCE_TYPE      = {'mrr','mtt_p_mpp','mtr','mtp'};
    HARD_DIR         = {'.MZZ','.MXX_P_MYY','.MXZ_MYZ','.MXY_MXX_M_MYY'}; % These are hard-coded directory names in AxiSEM
    NSOURCE          = 4;
    SOURCE_AMPLITUDE = repmat(event.Mo(:),1,4); % AxiSEM examples use this
else
    error(['Not setup for SIMULATION_TYPE = ',SIMULATION_TYPE]);
end

% Get name of 1D background model
theFile = [prjDir,'/AXISEM_INPUT/inparam_mesh'];
BACKGROUND_MODEL = sf_read_input(theFile,'BACKGROUND_MODEL','BACKGROUND_MODEL',false);
if strcmp(BACKGROUND_MODEL,'external')
    BACKGROUND_MODEL = sf_read_input(theFile,'EXT_MODEL','EXT_MODEL',false);
    warning('Did you fix AxiSEM to accept external models yet? If not, this simulation will error.')
else
    BACKGROUND_MODEL = 'iasp91_dsm';
end

% Derive geographic model dimensions
% Points defining cartesian lengths along box at center and ends
x0 = [xmax_m xmax_m xmax_m xmin_m 0 xmax_m]./1000;
y0 = [ymin_m 0 ymax_m ymax_m ymax_m ymax_m]./1000;
x1 = [xmin_m xmin_m xmin_m xmin_m 0 xmax_m]./1000;
y1 = [ymin_m 0 ymax_m ymin_m ymin_m ymin_m]./1000;
% Use surface widths...
z0      = zeros((size(x0)));
tf_elev = true; % In this case, z0 refers to elevation wrt Earth's surface
% % Or use bottom widths
% z0      = -depth_km*ones(size(x0));
% tf_elev = false; % In this case, z0 refers cartesian depth (negative) of box
% Map to geographic coordinates
[LON0,LAT0] = cart2geo(x0,y0,z0,lon0,lat0,tf_elev);
[LON1,LAT1] = cart2geo(x1,y1,z0,lon0,lat0,tf_elev);
% Calculate widths
dlon  = distance(LAT0,LON0,LAT1,LON1);
dlat  = max(dlon(4:6));
dlon  = max(dlon(1:3));
% Use maximum radial depth in model?
hr_km = depth_km;
% % Or minimum?
% tf_elev   = false; % In this case, z0 refers cartesian depth (negative) of box
% [~,~,ZR0] = cart2geo(x0,y0,-depth_km*ones(size(x0)),lon0,lat0,tf_elev);
% [~,~,ZR1] = cart2geo(x1,y1,-depth_km*ones(size(x0)),lon0,lat0,tf_elev);
% hr_km     = min([ZR0(:);ZR1(:)]-R_km);

%% Calculate minimum teleseismic P-wave travel-time to edge of box
% Used to guess at what time we start storing tractions along model edges
% Does this need to be the same across all simulations? See if reformat.par
% is required to run run_initialization.

% Find minimum distance from bottom edge of model to source (deg.)
xgeo             = [linspace(xmin_m,xmax_m,NX+1),xmax_m*ones(1,NY+1),linspace(xmin_m,xmax_m,NX+1),xmin_m*ones(1,NY+1)];
ygeo             = [ymax_m*ones(1,NX+1),linspace(ymin_m,ymax_m,NY+1),ymin_m*ones(1,NX+1),linspace(ymin_m,ymax_m,NY+1)];
tf_elev          = false; % In this case, z refers cartesian depth (negative) of box
[xgeo,ygeo,hgeo] = cart2geo(xgeo./1000,ygeo./1000,-depth_km*ones(size(xgeo)),lon0,lat0,tf_elev);
hgeo             = R_km - hgeo;

% Calculate minimum travel-times to model domain for each source
ttmin = zeros(nsrc,1);
tic
for ii = 1:nsrc
    d     = distance(event.latitude(ii),event.longitude(ii),ygeo,xgeo);
    imin  = find(d == min(d),1,'first');
    ttime = taup_time(a1Dmodel,aPhase_init,d(imin),event.depth(ii),hgeo(imin));
    % OLD SLOW TAUP SYSTEM CALL
    % ttime = call_tauP_time(a1Dmodel,num2str(d),'-h',num2str(event.depth(ii)),'--stadepth',num2str(hr_km),'-ph',aPhase_init);
    if isempty(ttime)
        error(['No arrival at ',num2str(d),' deg. for phase ',aPhase_init]);
    else
        ttmin(ii) = ttime;
    end
end
toc

% Use minimum travel-time to define time interval over which AxiSEM stores
% the tractions on the SPECFEM model
tmin = floor(ttmin) - dt_off;
tmax = tmin + DT*NSTEP;

% Calculating maximum travel-time to furthest station for each source. This
% step not strictly necessary but useful to check that phases of interest
% are captured in simulation.
disp('Calculating maximum travel-times for each source.');
ttmax = zeros(nsrc,1);
tic
for ii = 1:nsrc
    d     = max(distance(event.latitude(ii),event.longitude(ii),station.latitude,station.longitude));
    ttime = taup_time(a1Dmodel,aPhase_init,d,event.depth(ii),0);
    % OLD SLOW TAUP SYSTEM CALL
    % ttime = call_tauP_time(a1Dmodel,num2str(d),'-h',num2str(event.depth(ii)),'-ph',aPhase_fin);
    if isempty(ttime)
        warning(['No arrival at ',num2str(d),' deg. for phase ',aPhase_fin]);
    else
        ttmax(ii) = ttime;
    end
end
toc
if any((tmax-ttmax) < 0)
    warning('Simulation time interval may not capture all phases of interest.')
end

%% (1) WRITE PARFILEMESHCHUNK
theFile = [prjDir,'/MESH/ParFileMeshChunk'];

if isfile(theFile)
    error(['The file ',theFile,' already exists!']);
else
    fid = fopen(theFile,'w');
    
    fprintf(fid,'%s\n','# Domain Width: Longitude (deg.) Latitude (deg.)');
    fprintf(fid,'%9.5f %9.5f\n',dlon,dlat);
    fprintf(fid,'%s\n','# Domain Origin: Longitude (deg.) Latitude (deg.) Azimuth (deg.)');
    fprintf(fid,'%10.5f %10.5f %10.5f\n',lon0,lat0,azim);
    fprintf(fid,'%s\n','# Radial thickness: Depth (km)');
    fprintf(fid,'%10.5f\n',hr_km);
    fprintf(fid,'%s\n','# Number of elements: NX NY NZ');
    fprintf(fid,'%5.0f %5.0f %5.0f\n',NX,NY,NZ);
    fprintf(fid,'%s\n','# Background model: Name');
    fprintf(fid,'%s\n',BACKGROUND_MODEL);
    fprintf(fid,'%s\n','# Burried domain? Logical');
    fprintf(fid,'%s',subsur);
    
    fclose(fid);
end

%% (2) WRITE REFORMAT.PAR
theFile = [prjDir,'/source_files/reformat_'];
for ii = 1:nsrc
    if isfile([theFile,'CMTSOLUTION_',num2str(event.id(ii)),'.par'])
        error(['The file ',theFile,' already exists!']);
    else
        fid = fopen([theFile,'CMTSOLUTION_',num2str(event.id(ii)),'.par'],'w');
        
        fprintf(fid,'%f\n',1/DT);
        fprintf(fid,'%.0f %.0f',tmin(ii),tmax(ii));
        
        fclose(fid);
    end
end

%% (3) WRITE EXPAND_2D_3D
theFile = [prjDir,'/source_files/expand_2D_3D_'];

for ii = 1:nsrc
    apnd = ['CMTSOLUTION_',num2str(event.id(ii)),'.par'];
    
    if isfile([theFile,apnd])
        error(['The file ',[theFile,apnd],' already exists!']);
    end
        
    fid = fopen([theFile,apnd],'w');
    
    fprintf(fid,'%s\n','input_box.txt');
    fprintf(fid,'%s\n','input_box_sem_cart.txt');
    fprintf(fid,'%5.0f %s\n',NTHETA*NRADIAL,'  # Number of AxiSEM MPI processors used in solver');
    fprintf(fid,'%9.4f %9.4f %s\n',event.latitude(ii),event.longitude(ii),'  # Source location (Lat. Lon.)');
    fprintf(fid,'%10.5f %10.5f %10.5f %s\n',lat0,lon0,azim,'  # Geographic coordinates of model origin and orientation (Lat. Lon. azimuth)');
    fprintf(fid,'%s\n','.false.   # Used for KH integral and reciprocity; usually false');
    fprintf(fid,'%5.0f %s\n',NSOURCE,'  # Number of AxiSEM simulations (depends on source type)');
    fprintf(fid,'%5.0f %s\n',1,'  # A dummy value (just for KH)');
    fprintf(fid,'%5.0f %s\n',NPROC,'  # Number of SPECFEM MPI processors used in solver');
    % The relative paths of some directories relative to the run_axisem/SOLVER/RunName directory
    fprintf(fid,'%s\n','../../../MESH/');
    fprintf(fid,'%s\n',['../../../',LOCAL_PATH,'/']);
    fprintf(fid,'%s\n',['../../../',TRACTION_PATH,'/']);
    
    fclose(fid);
end

%% (4) WRITE INPARAM_SOURCE
theFile = [prjDir,'/source_files/inparam_source_'];

for jj = 1:nsrc
    for ii = 1:length(SOURCE_TYPE)
        apnd = cat(2,['CMTSOLUTION_',num2str(event.id(jj))],HARD_DIR{ii});
        
        if isfile([theFile,apnd])
            error(['The file ',[theFile,apnd],' already exists!']);
        else
            fid = fopen([theFile,apnd],'w');
            
            % Some information on source types
            fprintf(fid,'%s\n','# Basic source types');
            fprintf(fid,'%s\n','# For monopoles, specify one of the following: mrr, explosion, mtt_p_mpp, or vertforce');
            fprintf(fid,'%s\n','# For dipoles, specify one of the following: mtr, mpr, thetaforce, or phiforce');
            fprintf(fid,'%s\n','# For quadrupoles, specify one of the following: mtp or mtt_m_mpp');
            % Source parameters
            fprintf(fid,'%s\n',['SOURCE_TYPE       ',SOURCE_TYPE{ii}]);
            fprintf(fid,'\n%s\n','# Depth of source (km)');
            fprintf(fid,'%s %f\n','SOURCE_DEPTH       ',event.depth(jj));
            fprintf(fid,'\n%s\n','# Source latitude (deg.)');
            fprintf(fid,'%s %f\n','SOURCE_LAT       ',event.latitude(jj));
            fprintf(fid,'\n%s\n','# Source Longitude (deg.)');
            fprintf(fid,'%s %f\n','SOURCE_LON       ',event.longitude(jj));
            fprintf(fid,'\n%s\n','# Amplitude of source (Nm; N for force source)');
            fprintf(fid,'%s %E','SOURCE_AMPLITUDE       ',SOURCE_AMPLITUDE(jj,ii));
            
            fclose(fid);
        end
    end
end

%% (5) WRITE A DUMMY STATION FILE
% This is for AxiSEM...but I do not think it is really used for the coupled
% code. I think AxiSEM just needs it in the case that it needs to write out
% seismograms...but this is done by SPECFEM
%
% We could use the STATION file for SPECFEM but AxiSEM requires station
% locations in geographic coordinates and the coupled code works in
% cartesian so this is the simplest option for now.

theFile = [prjDir,'/AXISEM_INPUT/STATIONS_DUMMY'];

fid = fopen(theFile,'w');
fprintf(fid,'%5s %5s %10.5f %10.5f %6.1f %6.1f\n','101','DUM',lat0,lon0,0,0);
fclose(fid);

%% CHECK THESE VALUES FOR CONSISTANCY
% Can run this section alone if you define prjDir

% 1) Mesh_Par_file
% --> NEX_XI and NEX_ETA must be divisible by NPROC_XI and NPROC_ETA,
%     respectively
theFile   = [prjDir,'/DATA/Par_file'];
NPROC     = sf_read_input(theFile,'NPROC','=',true);
theFile   = [prjDir,'/DATA/meshfem3D_files/Mesh_Par_file'];
NEX_XI    = sf_read_input(theFile,'NEX_XI','=',true);
NEX_ETA   = sf_read_input(theFile,'NEX_ETA','=',true);
NPROC_XI  = sf_read_input(theFile,'NPROC_XI','=',true);
NPROC_ETA = sf_read_input(theFile,'NPROC_ETA','=',true);
if ~(rem(NEX_XI,NPROC_XI) == 0)
    warning('NEX_XI is not divisible by NPROC_XI!');
end
if ~(rem(NEX_ETA,NPROC_ETA) == 0)
    warning('NEX_ETA is not divisible by NPROC_ETA!');
end
% --> NPROC_XI * NPROC_ETA must equal NPROC in Par_file
if ~(NPROC_XI*NPROC_ETA == NPROC) && ~(NPROC_XI*NPROC_ETA == 1)
    warning('NPROC_XI*NPROC_ETA does NOT equal NPROC')
end

% % 2) Possibly problematic if tractions from AxiSEM end before SPECFEM
% % simulation (depends on phases of interest)
% theFile           = [prjDir,'/AXISEM_INPUT/inparam_basic'];
% SEISMOGRAM_LENGTH = sf_read_input(theFile,'SEISMOGRAM_LENGTH','SEISMOGRAM_LENGTH',true);
% tmax              = dlmread([prjDir,'/AXISEM_INPUT/reformat.par']);
% tmax              = tmax(end);
% if tmax > SEISMOGRAM_LENGTH
%     warning('SPECFEM simulation time exceeds AxiSEM simulation time. May need to increase SEISMOGRAM_LENGTH in inparam_basic.');
% else
%     disp('VA BENE!');
% end

% 3) expand_2D_3D.par
% --> Number of processors correctly defined for AxiSEM/SPECFEM
% --> Should be okay if generated with this script
