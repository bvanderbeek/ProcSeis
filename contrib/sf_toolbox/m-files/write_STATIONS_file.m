%% Create Synthetic Station File
% Creates a STATIONS file for an array of receiver locations
% Elevation Note:
% I think SPECFEM only uses the receiver burial depths NOT the elevations.
% The burial depth gives the position of the receiver with respect to the
% surface of the mesh. If USE_SOURCES_RECEIVERS_Z is true, the burial depth
% is taken as the true cartesian position in the mesh coordinate system.
% For my purposes, keep USE_SOURCES_RECEIVERS_Z false.
%
% All stations assumed to be on surface in this script
clear; close all; clc

% INPUT
theFile   = '/Users/bvanderbeek/research/NEWTON/AxiSEM_SPECFEM3D_Projects/examples/local_source/DATA/STATIONS'; % File that will be written
slat      = (-100:50:100)*1000; % Latitude (dd) or y-coordinate (m)
slon      = (-100:50:100)*1000; % Longitude (dd) or x-coordinate (m)
NW        = 'SYN'; % A network code

%% Define station parameters

% Derived station parameters
nsta      = length(slat)*length(slon);
ID        = (10^ceil(log10(nsta)))+linspace(1,nsta,nsta); % Station id/name. MUST BE INTEGER for this script
[LAT,LON] = meshgrid(slat,slon);
Z         = zeros(size(LAT)); % Elevation
D         = zeros(size(LAT)); % Depth
% Vectorize
LAT = LAT(:);
LON = LON(:);
Z   = Z(:);
D   = D(:);
ID  = ID(:);

%% Write file

% Create file
if isfile(theFile)
    error(['The file, ',theFile,' already exists!'])
else
    fid = fopen(theFile,'w');
end
for ii = 1:length(ID)
    fprintf(fid,'%6s %6s %14.5f %14.5f %10.3f %10.3f\n',num2str(ID(ii)),NW,LAT(ii),LON(ii),Z(ii),D(ii));
end
fclose(fid);
