%% Write multiple source files
% Creates multiple CMTSOLUTION files for a random distribution of sources.
clear; close all; clc

% INPUT
outDir  = '/Users/bvanderbeek/research/NEWTON/AxiSEM_SPECFEM3D_Projects/examples/tele_double_couple/source_files';
lat0    = 0; % Origin of array latitude (dd)
lon0    = 0; % Origin of array longitude (dd)
Mw      = 6.5; % Magnitude
t       = 10; % Half-duration (s)
strike  = 0; % Fault strike relative to ray azimuth from source origin
dip     = 60; % Fault dip (deg.)
rake    = 90; % Fault rake (deg.)
frot    = 10; % Rotation of fault plane with respect to source-receiver azimuth (+CW from North; deg.)
% Parameters for a random distribution of sources
N       = 8; % Number of event files to create
dlimits = [40 90]; % Min and max range (deg.)
alimits = [0 360]; % Min and max azimuth (deg.)
hlimits = [10 100]; % Min and max depth (km)

%% Create Random Distribution of sources

% Generate a random arrays of source distances, back-azimuths, and depths
D   = dlimits(1) + (dlimits(2)-dlimits(1))*rand(N,1);
baz = alimits(1) + (alimits(2)-alimits(1))*rand(N,1);
H   = hlimits(1) + (hlimits(2)-hlimits(1))*rand(N,1);
% Synthetic source coordinates
[slat,slon] = reckon(lat0,lon0,D,baz);
% Source-to-receiver azimuth
[~,az]  = distance(slat,slon,lat0,lon0);

% Plot
polarplot(baz*pi./180,D,'.b','markersize',15);
title('Source Back-Azimuths and Ranges');

%% Write the files

for ii = 1:N
    % Output file
    theFile = [outDir,'/CMTSOLUTION_',num2str(ii)];
    
    % Rotate fault and write file
    gamma = az(ii) - 180 + frot; % Azimuth of fault plane + CW from N.
    sf_write_CMTSOLUTION_file(theFile,slon(ii),slat(ii),H(ii),t,Mw,[],gamma+strike,dip,rake);
end
