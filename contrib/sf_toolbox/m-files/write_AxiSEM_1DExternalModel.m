%% Write AxiSEM External 1D Model file
% This creates a simple 1D velocity model file that can be used with
% AxiSEM. This is useful if you do not want to use the built-in 1D models.
%
% I am a bit confused about this file type. It seems to store coefficients
% that describe velocity variations within each layer but I am not sure of
% the form of the equation that these coefficients describe. For now, I
% just create many layers and assign their velocity to the first
% coefficient and leave the remaining zero.
%
% Can do anelastic and anisotropic models but hard-coded now for elastic
% isotropic models.
%
% Note that the AxiSEM release included with SPECFEM3D is hard-coded for
% models with 12 layers. To use an arbitrary number of layers the AxiSEM
% mesh routines must be modified and re-compiled. This is easy but
% annoying.
clear; close all; clc

% INPUT
theFile = 'RadialEarthModels/dsm131_custom_tomography_model_karatowu2_10km_resol_td_isotropic_fossil_i1j1'; % Output file name
theData = 'RadialEarthModels/iasp91_R_rho_vp_vs.mat'; % Radial Earth model MAT-file with R, rho, vp, and vs variables (km, km/s); discontinuities should be built in via repeated depth nodes
theRef  = 'RadialEarthModels/data_tomography_model_karatowu2_10km_resol_td_isotropic_fossil_refi1j1.mat'; % If empty, then just uses IASP91
tf_NoCrust = false; % If true, removes crustal velocities

%% Write model file

% Create file
if exist(theFile,'file') > 0
    error(['The file ', theFile, 'already exists.']);
else
    fid = fopen(theFile,'w');
end

% Load model data
load(theData);

figure; hold on;
plot(vp,R,'-b');
plot(vs,R,'-r');
plot(rho,R,'-g');
box on; grid on;

% Remove duplicate points
idup        = find(diff(R) == 0);
rho(idup)   = (rho(idup) + rho(idup+1))./2;
vp(idup)    = (vp(idup) + vp(idup+1))./2;
vs(idup)    = (vs(idup) + vs(idup+1))./2;
R(idup+1)   = [];
rho(idup+1) = [];
vp(idup+1)  = [];
vs(idup+1)  = [];

% Any radial earth model may now be used but it must contain 131 points
% exactly (or interpolate to the radial depth defined above).
if ~isempty(theRef)
    % Load reference model structure and interpolate
    load(theRef);
    data.vp  = interp1(data.r,data.vp,R);
    data.vs  = interp1(data.r,data.vs,R);
    data.rho = interp1(data.r,data.rho,R);
    
    % Reset velocity profile values
    reset      = ~isnan(data.vp);
    vp(reset)  = data.vp(reset);
    vs(reset)  = data.vs(reset);
    rho(reset) = data.rho(reset);
    
    plot(vp,R,'--k');
    plot(vs,R,'--k');
    plot(rho,R,'--m');
    shg;
end

% Reset crustal velocities
if tf_NoCrust
    imoho = find(vp <= 6.5,1,'first');
    vp(imoho:end)  = vp(imoho-1);
    vs(imoho:end)  = vs(imoho-1);
    rho(imoho:end) = rho(imoho-1);
end

% Write the header
nzone = length(R);
fprintf(fid,'%3.0f  %s',nzone-1,'nzone');

% Write model parameters
nline = 0;
for ii = 1:(nzone-1)
    % There are 6 lines per layer
    % Line 1: R_min R_max rho ? ? ?
    val = (rho(ii) + rho(ii+1))/2;
    fprintf(fid,'\n%6.1f  %6.1f  %8.5f%11.5f%11.5f%11.5f',R(ii),R(ii+1),val,0,0,0);
    % Line 2: vpv ? ? ?
    val = (vp(ii) + vp(ii+1))/2;
    fprintf(fid,'\n%24.5f%11.5f%11.5f%11.5f',val,0,0,0);
    % Line 3: vph ? ? ?
    fprintf(fid,'\n%24.5f%11.5f%11.5f%11.5f',val,0,0,0);
    % Line 4: vsv ? ? ?
    val = (vs(ii) + vs(ii+1))/2;
    fprintf(fid,'\n%24.5f%11.5f%11.5f%11.5f',val,0,0,0);
    % Line 5: vsh ? ? ?
    fprintf(fid,'\n%24.5f%11.5f%11.5f%11.5f',val,0,0,0);
    % Line 6: eta ? ? ? qm qk
    fprintf(fid,'\n%24.5f%11.5f%11.5f%11.5f%6s%5s',1,0,0,0,'0.d0','0.d0');
    
    nline = nline + 1;
end

fclose(fid);
