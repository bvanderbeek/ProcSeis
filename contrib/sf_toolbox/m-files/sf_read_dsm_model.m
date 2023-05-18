function [R,rho,Vph,Vpv,Vsh,Vsv,eta] = sf_read_dsm_model(theFile)

% Open file and read header
fid   = fopen(theFile);
nline = strsplit(fgetl(fid));
nmax  = str2double(nline{1}); % Number of layers
n     = 1; % Initialize layer counter

% Loop over layers
R   = zeros(nmax,2);
rho = zeros(nmax,4);
Vph = zeros(nmax,4);
Vpv = zeros(nmax,4);
Vsh = zeros(nmax,4);
Vsv = zeros(nmax,4);
eta = zeros(nmax,4);
while n <= nmax
    % Read first line of layer
    nline = strsplit(strtrim(fgetl(fid)));
    % Top and bottom of layer
    R(n,1) = str2double(nline{1});
    R(n,2) = str2double(nline{2});
    % Density coefficients
    rho(n,1) = str2double(nline{3});
    rho(n,2) = str2double(nline{4});
    rho(n,3) = str2double(nline{5});
    rho(n,4) = str2double(nline{6});
    
    % Read second line of layer
    nline = strsplit(strtrim(fgetl(fid)));
    % Vpv coefficients
    Vpv(n,1) = str2double(nline{1});
    Vpv(n,2) = str2double(nline{2});
    Vpv(n,3) = str2double(nline{3});
    Vpv(n,4) = str2double(nline{4});
    
    % Read third line of layer
    nline = strsplit(strtrim(fgetl(fid)));
    % Vpv coefficients
    Vph(n,1) = str2double(nline{1});
    Vph(n,2) = str2double(nline{2});
    Vph(n,3) = str2double(nline{3});
    Vph(n,4) = str2double(nline{4});
    
    % Read fourth line of layer
    nline = strsplit(strtrim(fgetl(fid)));
    % Vpv coefficients
    Vsv(n,1) = str2double(nline{1});
    Vsv(n,2) = str2double(nline{2});
    Vsv(n,3) = str2double(nline{3});
    Vsv(n,4) = str2double(nline{4});
    
    % Read fith line of layer
    nline = strsplit(strtrim(fgetl(fid)));
    % Vpv coefficients
    Vsh(n,1) = str2double(nline{1});
    Vsh(n,2) = str2double(nline{2});
    Vsh(n,3) = str2double(nline{3});
    Vsh(n,4) = str2double(nline{4});
    
    % Read sixth line of layer
    nline = strsplit(strtrim(fgetl(fid)));
    % Vpv coefficients
    eta(n,1) = str2double(nline{1});
    eta(n,2) = str2double(nline{2});
    eta(n,3) = str2double(nline{3});
    eta(n,4) = str2double(nline{4});
    % Shear and bulk attenuation...funny format; read in as NaN
    % qm     = str2double(nline{5});
    % qk     = str2double(nline{6});
    
    % Update counters
    n = n + 1;
end
fclose(fid);
