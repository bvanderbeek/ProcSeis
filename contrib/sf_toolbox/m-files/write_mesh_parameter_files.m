%% Create Mesh Parameter Files
% This script generates files that are required for creating meshes with
% SPECFEM's internal mesher. This script will also create the required
% simulation directories in the project folder. When setting up a new
% project, this should probably be the first script you run.
%
% Files created:
% 1) Mesh_Par_file
% 2) interfaces.dat
% 3) interfaceX (the z-coordinates of planar interfaces defined in interfaces.dat
%
% Notes:
% + Parameters relevant to C-PML absorbing boundary conditions are
%   currently hard-coded. This is because the internal mesher cannot use
%   these types of boundary conditions and so these values are not used.
%
% B. VanderBeek (6-April-2019)
clear; close all; clc

% Add sf_toolbox
addpath(genpath('/Users/bvanderbeek/research/software/sf_toolbox/m-files'));
% Add path to TauP (uses IASP91 velocity model to estimate resolution)
setenv('TAUPJAR','/Users/bvanderbeek/research/software/sf_toolbox/TauP/TauP-2.4.5.jar');
addpath(genpath('/Users/bvanderbeek/research/software/sf_toolbox/TauP/m-files'));

% Location of simulation directory (i.e. from where SPECFEM is run)
prjDir   = ''; %'/Users/bvanderbeek/research/NEWTON/AxiSEM_SPECFEM3D_Projects/projects/SubductionModel_101_Centered_Annulus_ISO';

% Lateral mesh dimensions
xmin     = -1000000; % Minimum longitude (dd) or x-coordinate (m)
xmax     = 1000000; % Maximum longitude (dd) or x-coordinate (m)
ymin     = -1500000; % Minimum latitude (dd) or y-coordinate (m)
ymax     = 1500000; % Maximum latitude (dd) or y-coordinate (m)
% Ideally you want elements to be square to reduce errors in the solution.
% However, SPECFEM requires that the number of elements be divisible by the
% number of proccessor used in the simulation. These values are just
% initial estimates and will be adjusted as needed by this script.
nx       = 128; % Desired number of *elements* in x-direction (--)
ny       = 192; % Desired number of *elements* in y-direction (--)
npx      = 8; % Number of MPI chunks in the x-direction (--)
npy      = 8; % Number of MPI chunks in the y-direction (--)

% Vertical mesh dimension
% For SPECFEM, depth is negative into the earth and zmin is the
% z-coordinate of the mesh bottom in the Earth. In this script, The surface
% of the mesh is assumed to coincide with the Earth's surface (z = 0) but
% this does not need to be the case.
% If nz is left empty, this script will chose the total number of elements 
% in the z-direction such that elements remain as square as possible.
zmin     = -1000000; % Minimum z-coordinate (m); depth is *negative* into earth!
nz       = []; % Desired number of elements in z-direction (--)
tf_Curve = true; % If true, correct cartesian grid and interfaces for Earth's curvature (logical)

% Planar mesh interfaces. Do *not* include free surface (z = 0);
% Define interface z-coordinates sequentially with increasing depth.
% Interface depths (including the mesh surface) will be corrected for 
% Earth's curvature if tf_Curv is true. If doubling across these layers, it
% is useful to put the interafce a little below the true interface
zint     = -120000; % z-coordinate of interface (m); depth is *negative* into earth!
tf_dbl   = true; % If true, the element size will double across these interfaces (logical)

% UTM parameters
% Unless mesh dimensions are given in UTM coordinates, tf_NoUTM should be
% true. Note that tf_NoUTM must be true for use with AxiSEM.
tf_NoUTM = true; % If true, surpresses UTM project (logical). Should be true for AxiSEM coupling.
UTM      = []; % The UTM projection zone (--); no value need if tf_NoUTM is true

% If SAVE_MESH_FILES is true in the Par_file, then mesh files are output in
% these formats where true
tf_abq = false; % (logical)
tf_dxf = false; % (logical)
tf_vtk = false; % (logical)


%% Define z elements
% Hard-coded location of DATABASES_MPI directory
DBdir = './DATABASES_MPI';

% Force number of elements in x,y to be divisible by number of processors
if tf_dbl
    nx = 8*npx*round(nx/(8*npx));
    ny = 8*npy*round(ny/(8*npy));
else
    nx = npx*round(nx/npx);
    ny = npy*round(ny/npy);
end

% Average lateral grid-spacing (ideally equal)
d  = (((xmax-xmin)/nx) + ((ymax-ymin)/ny))/2;

% Number of potential doubling layers
ndbl = length(zint);
zint = cat(1,0,zint(:),zmin);
zint = flipud(zint); % Interfaces defined starting with bottom

% Get number of elements within layers
if (ndbl > 0) && tf_dbl
    % Construct system of equations
    A    = eye(ndbl+1);
    eqnz = zeros(1,ndbl+1);
    for ii = 1:ndbl
        A(ii,ndbl+1) = -2^(ii-ndbl-1);
        eqnz(ii)     = zint(ii+1)-zint(ii);
    end
    eqnz(ii+1) = zint(ii+2) - zint(ii+1);
    % Add constraint equation
    if isempty(nz)
        % Square-element constraint
        b = [zeros(ndbl,1); 1/d];
    else
        % Number of elements in z-direction constraint
        A(end,:) = eqnz;
        b = [zeros(ndbl,1); nz];
    end
    % Solve system
    idz = lsqr(A,b,[],1000);
    % Get number of nodes in z-direction within each layer
    N   = round(eqnz(:).*idz);
    nz  = sum(N);
    dz  = eqnz(:)./N(:);
else
    if isempty(nz)
        N  = round(-zmin/d);
        nz = N;
    else
        N = nz;
    end
    dz = -zmin/N;
end

%% Create interfaces
R       = 6371*1000;
xg      = linspace(xmin,xmax,nx+1);
yg      = linspace(ymin,ymax,ny+1);
[Yi,Xi] = meshgrid(yg,xg);
Zi      = zeros(nx+1,ny+1,ndbl+1);

% Interface grid spacing
dx = (xmax-xmin)/nx;
dy = (ymax-ymin)/ny;

for ii = 1:(ndbl+1)
    % Index is ii + 1 to skip bottom of model
    % Correct for Earth curvature
    if tf_Curve
        Zi(:,:,ii) = sqrt(((R+zint(ii+1))^2) - (Xi.^2) - (Yi.^2)) - R;
    else
        Zi(:,:,ii) = zint(ii+1)*ones(nx+1,ny+1);
    end
end

disp('Final Model Dimensions at Surface:');
disp(['NX = ',num2str(nx),'; NY = ',num2str(ny),'; NZ = ',num2str(nz)]);
disp(['DX = ',num2str(dx/1000),'; DY = ',num2str(dy/1000),'; DZ = ',num2str(dz(end)/1000),' (km)']);

%% Make required directories

system(['mkdir -p ',prjDir,'/DATA/meshfem3D_files']);

%% Write Mesh_Par_file

if isfile([prjDir,'/DATA/meshfem3D_files/Mesh_Par_file'])
    error(['The file ',[prjDir,'/DATA/meshfem3D_files/Mesh_Par_file'],' already exists!']);
else
    fid = fopen([prjDir,'/DATA/meshfem3D_files/Mesh_Par_file'],'w');
end

fprintf(fid,'%s\n','#-----------------------------------------------------------');
fprintf(fid,'%s\n','#');
fprintf(fid,'%s\n','# M E S H   I N P U T   P A R A M E T E R S');
fprintf(fid,'%s\n','#');
fprintf(fid,'%s\n','#-----------------------------------------------------------');

fprintf(fid,'\n%s\n','# Dimensions of mesh block');
fprintf(fid,'%s\n',['LATITUDE_MIN                    = ',num2str(ymin,'%.2f')]);
fprintf(fid,'%s\n',['LATITUDE_MAX                    = ',num2str(ymax,'%.2f')]);
fprintf(fid,'%s\n',['LONGITUDE_MIN                   = ',num2str(xmin,'%.2f')]);
fprintf(fid,'%s\n',['LONGITUDE_MAX                   = ',num2str(xmax,'%.2f')]);
fprintf(fid,'%s\n',['DEPTH_BLOCK_KM                  = ',num2str(abs(zmin)/1000,'%.2f')]); % The only paramter in kilometers
if tf_NoUTM
    fprintf(fid,'%s\n','UTM_PROJECTION_ZONE             = 11'); % Not used if tf_NoUTM is true
    fprintf(fid,'%s\n','SUPPRESS_UTM_PROJECTION         = .true.');
else
    fprintf(fid,'%s\n',['UTM_PROJECTION_ZONE             = ',num2str(UTM)]);
    fprintf(fid,'%s\n','SUPPRESS_UTM_PROJECTION         = .false.');
end

fprintf(fid,'\n%s\n','# Name of interface file');
fprintf(fid,'%s\n','INTERFACES_FILE                 = interfaces.dat');

fprintf(fid,'\n%s\n','# Name of cavity file');
fprintf(fid,'%s\n','CAVITY_FILE                     = NoCavity.dat'); % DUMMY VALUES

fprintf(fid,'\n%s\n','# Number of elements along x(xi)- and y(eta)-directions at surface');
fprintf(fid,'%s\n','# Must be a multiple of NPROC_* below');
fprintf(fid,'%s\n',['NEX_XI                          = ',num2str(nx)]);
fprintf(fid,'%s\n',['NEX_ETA                         = ',num2str(ny)]);

fprintf(fid,'\n%s\n','# Number of MPI processes along x(xi)- and y(eta)-directions');
fprintf(fid,'%s\n',['NPROC_XI                        = ',num2str(npx)]);
fprintf(fid,'%s\n',['NPROC_ETA                       = ',num2str(npy)]);


fprintf(fid,'\n%s\n','#-----------------------------------------------------------');
fprintf(fid,'%s\n','#');
fprintf(fid,'%s\n','# D O U B L I N G   L A Y E R S');
fprintf(fid,'%s\n','#');
fprintf(fid,'%s\n','#-----------------------------------------------------------');

fprintf(fid,'\n%s\n','# Use irregular mesh (i.e. with doublings)?');
if tf_dbl
    fprintf(fid,'%s\n',['USE_REGULAR_MESH                = ','.false.']);
else
    fprintf(fid,'%s\n',['USE_REGULAR_MESH                = ','.true.']);
end

if (ndbl > 0) && tf_dbl
    fprintf(fid,'\n%s\n','# Number of doubling layers');
    fprintf(fid,'%s\n',['NDOUBLINGS                      = ',num2str(ndbl)]);
    
    fprintf(fid,'\n%s\n','# Index of doubling layers (counts from bottom up)');
    idouble = 0;
    for ii = 1:ndbl
        idouble = idouble + N(ii);
        fprintf(fid,'%s\n',['NZ_DOUBLING_',num2str(ii),'                   = ',num2str(idouble)]);
    end
else
    fprintf(fid,'\n%s\n','# Number of doubling layers');
    fprintf(fid,'%s\n','NDOUBLINGS                      = 1');
    fprintf(fid,'\n%s\n','# Index of doubling layers (counts from bottom up)');
    fprintf(fid,'%s\n','NZ_DOUBLING_1                   = 1');
end


fprintf(fid,'\n%s\n','#-----------------------------------------------------------');
fprintf(fid,'%s\n','#');
fprintf(fid,'%s\n','# V I S U A L I Z A T I O N');
fprintf(fid,'%s\n','#');
fprintf(fid,'%s\n','#-----------------------------------------------------------');

fprintf(fid,'\n%s\n','# Export mesh files in following formats?');
if tf_abq
    fprintf(fid,'%s\n','CREATE_ABAQUS_FILES             = .true.');
else
    fprintf(fid,'%s\n','CREATE_ABAQUS_FILES             = .false.');
end
if tf_dxf
    fprintf(fid,'%s\n','CREATE_DX_FILES                 = .true.');
else
    fprintf(fid,'%s\n','CREATE_DX_FILES                 = .false.');
end
if tf_vtk
    fprintf(fid,'%s\n','CREATE_VTK_FILES                = .true.');
else
    fprintf(fid,'%s\n','CREATE_VTK_FILES                = .false.');
end

fprintf(fid,'\n%s\n','# Path to store database files');
fprintf(fid,'%s\n',['LOCAL_PATH                      = ',DBdir]);


fprintf(fid,'\n%s\n','#-----------------------------------------------------------');
fprintf(fid,'%s\n','#');
fprintf(fid,'%s\n','# C M P L');
fprintf(fid,'%s\n','#');
fprintf(fid,'%s\n','#-----------------------------------------------------------');

fprintf(fid,'\n%s\n','# CPML perfectly matched absorbing layers');
fprintf(fid,'%s\n','THICKNESS_OF_X_PML              = 12.3d0');
fprintf(fid,'%s\n','THICKNESS_OF_Y_PML              = 12.3d0');
fprintf(fid,'%s\n','THICKNESS_OF_Z_PML              = 12.3d0');


fprintf(fid,'\n%s\n','#-----------------------------------------------------------');
fprintf(fid,'%s\n','#');
fprintf(fid,'%s\n','# D O M A I N   M A T E R I A L S');
fprintf(fid,'%s\n','#');
fprintf(fid,'%s\n','#-----------------------------------------------------------');

fprintf(fid,'\n%s\n','# Number of materials');
fprintf(fid,'%s\n','NMATERIALS                      = 1'); % DUMMY VALUES
fprintf(fid,'%s\n','# Definition of material properties is as follows:');
fprintf(fid,'%s\n','# #material_id  #rho  #vp  #vs  #Q_Kappa  #Q_mu  #anisotropy_flag  #domain_id');
fprintf(fid,'%s\n','#     Q_Kappa          : Q_Kappa attenuation quality factor');
fprintf(fid,'%s\n','#     Q_mu             : Q_mu attenuation quality factor');
fprintf(fid,'%s\n','#     anisotropy_flag  : 0 = no anisotropy / 1,2,... check the implementation in file aniso_model.f90');
fprintf(fid,'%s\n','#     domain_id        : 1 = acoustic / 2 = elastic');
fprintf(fid,'%s\n','1 2000 3500 2000 9999. 9999. 1 2'); % DUMMY VALUES


fprintf(fid,'\n%s\n','#-----------------------------------------------------------');
fprintf(fid,'%s\n','#');
fprintf(fid,'%s\n','# D O M A I N   R E G I O N S');
fprintf(fid,'%s\n','#');
fprintf(fid,'%s\n','#-----------------------------------------------------------');

fprintf(fid,'\n%s\n','# Number of regions');
fprintf(fid,'%s\n','NREGIONS                        = 1'); % DUMMY VALUES
fprintf(fid,'%s\n','# Define each model region as follows:');
fprintf(fid,'%s\n','# #NEX_XI_BEGIN  #NEX_XI_END  #NEX_ETA_BEGIN  #NEX_ETA_END  #NZ_BEGIN #NZ_END  #material_id');
fprintf(fid,'%.0f %.0f %.0f %.0f %.0f %.0f %.0f\n',1,nx,1,ny,1,nz,1); % DUMMY VALUES

fclose(fid);

%% Write the interface file

if isfile([prjDir,'/DATA/meshfem3D_files/interfaces.dat'])
    error(['The file ',[prjDir,'/DATA/meshfem3D_files/interfaces.dat'],' already exists!']);
else
    fid = fopen([prjDir,'/DATA/meshfem3D_files/interfaces.dat'],'w');
end

fprintf(fid,'%s\n','# Number of interfaces');
fprintf(fid,'%s\n',num2str(ndbl+1));

fprintf(fid,'%s\n','# Define parameters for each interface starting from the bottom as follows:');
fprintf(fid,'%s\n','# SUPPRESS_UTM_PROJECTION  NXI  NETA LONG_MIN   LAT_MIN    SPACING_XI SPACING_ETA');
fprintf(fid,'%s\n','# name_of_interface_file');
if tf_NoUTM
    mystr = '.true.';
else
    mystr = '.false.';
end
for ii = 1:(ndbl+1)
    fprintf(fid,'%s\n',['# Interface ',num2str(ii)]);
    fprintf(fid,'%s\n',[mystr,' ',num2str(nx+1),' ',num2str(ny+1),' ',num2str(xmin,'%.2f'),' ',num2str(ymin,'%.2f'),' ',num2str(dx,'%.2f'),' ',num2str(dy,'%.2f')]);
    fprintf(fid,'%s\n',['interface',num2str(ii)]);
    
    % Also, write the interface file too
    if isfile([prjDir,'/DATA/meshfem3D_files/interface',num2str(ii)])
        error(['The file ',[prjDir,'/DATA/meshfem3D_files/interface',num2str(ii)],' already exists!']);
    else
        data = Zi(:,:,ii);
        dlmwrite([prjDir,'/DATA/meshfem3D_files/interface',num2str(ii)],data(:),'precision','%11.3f');
    end
end
fprintf(fid,'%s\n','#');

fprintf(fid,'%s\n','# Define number of elements in vertical direction within each layer starting from bottom');
for ii = 1:(ndbl+1)
    fprintf(fid,'%s\n',['# Layer ',num2str(ii)]);
    fprintf(fid,'%s\n',num2str(N(ii)));
end

fclose(fid);

%% Check resolution

% Model node z-coordinate vector
zq = zmin + dz(1)*linspace(0,N(1),N(1)+1)';
for ii = 2:length(N)
    zq = cat(1,zq,max(zq) + dz(ii)*linspace(1,N(ii),N(ii))');
end
zq = flipud(zq)./1000;

% Model element dimensions
dx_vec = dx*ones(N(1)+1,1);
dy_vec = dy*ones(N(1)+1,1);
dz_vec = dz(1)*ones(N(1)+1,1);
for ii = 2:length(N)
    dx_vec = cat(1,dx_vec,ii*dx*ones(N(ii),1));
    dy_vec = cat(1,dy_vec,ii*dy*ones(N(ii),1));
    dz_vec = cat(1,dz_vec,dz(ii)*ones(N(ii),1));
end
dx_vec = flipud(dx_vec)./1000;
dy_vec = flipud(dy_vec)./1000;
dz_vec = flipud(dz_vec)./1000;

% The 1D velocity model
theFile = [prjDir,'/MESH/iasp91_dsm'];
if isfile(theFile)
    [R,~,vp,~,vs,~,~] = sf_read_dsm_model(theFile);
    R  = flipud(mean(R,2) - 6371);
    vp = flipud(vp(:,1));
    vs = flipud(vs(:,1));
    % Interpolate to model z-coordinate
    vp = interp1(R,vp,zq,'linear',vp(1));
    vs = interp1(R,vs,zq,'linear',vs(1));
else
    [vp,vs] = get_TauP_model(abs(zq),'iasp91');
end

% The approximate minimum period resolved for P and S waves (see manual).
% This is only estimate. Checker mesher output logs for more accurate
% estimates.
ngll = 5; % Number of gll integration points in each direction. This is fixed in SPECFEM when compiling (should be 5 unless changed)
Tp   = max([ngll*dx_vec./(vp*(ngll - 1)), ngll*dy_vec./(vp*(ngll - 1)), ngll*dz_vec./(vp*(ngll - 1))],[],2);
Ts   = max([ngll*dx_vec./(vs*(ngll - 1)), ngll*dy_vec./(vs*(ngll - 1)), ngll*dz_vec./(vs*(ngll - 1))],[],2);
disp(['Minimum period resolved using compressional velocity is ~',num2str(max(Tp)),' s.']);
disp(['Minimum period resolved using shear velocity is ~',num2str(max(Ts)),' s.']);

% Maximum time step allowed
dim     = 3; % Model dimensions
dt_crit = sf_critical_timestep(repmat(vp(:),3,1),[dx_vec(:);dy_vec(:);dz_vec(:)],ngll,dim);
dt_crit = min(dt_crit);
disp(['Estimated critical time-step is ',num2str(dt_crit),' s.']);
