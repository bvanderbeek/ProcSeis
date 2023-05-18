function [rho,vp,vs] = sf_get_axisem_1D(x,y,z,theFile)
% This is a matlab implementation of the routines found in
% generate_databases/model_coupled.f90 for defining the 1D model for
% coupling with AxiSEM.
%   + Direct translation; could be more efficient

% Earth radius
Re = 6371;

% Get default model file (in sf_toolbox) if none provided
if isempty(theFile)
    toolDir = which('sf_read_dsm_model');
    toolDir = strsplit(toolDir,'m-files');
    theFile = [toolDir{1},'RadialEarthModels/iasp91_dsm_original_11layer'];
end

% Read the AxiSEM model file
[zlayer,rho1D,vp1D,~,vs1D] = sf_read_dsm_model(theFile);

% Loop over interpolation points
nlayer = size(zlayer,1);
num    = length(x);
rho    = zeros(num,1);
vp     = zeros(num,1);
vs     = zeros(num,1);
for n = 1:num
    % Find layer index
    r = sqrt(x(n).^2 + y(n).^2 + (z(n)+Re).^2);
    il = 1;
    while (r > zlayer(il,1)) && (il < nlayer)
        il = il + 1;
    end
    il = il - 1;
    
    % Get model values
    r      = r/zlayer(nlayer,1);
    rho(n) = rho1D(il,1) + rho1D(il,2)*r + rho1D(il,3)*(r^2) + rho1D(il,4)*(r^3);
    vp(n)  = vp1D(il,1) + vp1D(il,2)*r + vp1D(il,3)*(r^2) + vp1D(il,4)*(r^3);
    vs(n)  = vs1D(il,1) + vs1D(il,2)*r + vs1D(il,3)*(r^2) + vs1D(il,4)*(r^3);
end
