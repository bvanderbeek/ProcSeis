function [H,data] = sf_read_ascii_model(theFile,varargin)

% Tomography model file is assumed to be on regular grid

% For fully anisotropic model this takes about 10 s per 1 million lines.

if isempty(varargin)
    tf_flip = false;
    tf_flat = false;
    ireset  = [];
elseif (length(varargin) == 1)
    tf_flip = varargin{1};
    tf_flat = false;
    ireset  = [];
elseif (length(varargin) == 2)
    tf_flip = varargin{1};
    tf_flat = varargin{2};
    ireset  = [];
elseif (length(varargin) == 3)
    tf_flip = varargin{1};
    tf_flat = varargin{2};
    ireset  = varargin{3};
else
    error('Problematic variable input arguments.')
end

% Open the file
fid = fopen(theFile);

% First line of header: model boundaries
char1 = '#';
nhead = -1; % track total number of header lines
while strcmp(char1,'#')
    nline = fgetl(fid);
    nline = strtrim(nline);
    char1 = nline(1);
    nhead = nhead + 1;
end
nline = strsplit(nline);
H.xmin = str2double(nline{1});
H.ymin = str2double(nline{2});
H.zmin = str2double(nline{3});
H.xmax = str2double(nline{4});
H.ymax = str2double(nline{5});
H.zmax = str2double(nline{6});

% Second line of header: model spacing
char1 = '#';
while strcmp(char1,'#')
    nline = fgetl(fid);
    nline = strtrim(nline);
    char1 = nline(1);
    nhead = nhead + 1;
end
nline = strsplit(nline);
H.dx = str2double(nline{1});
H.dy = str2double(nline{2});
H.dz = str2double(nline{3});

% Third line of header: model dimensions
char1 = '#';
while strcmp(char1,'#')
    nline = fgetl(fid);
    nline = strtrim(nline);
    char1 = nline(1);
    nhead = nhead + 1;
end
nline = strsplit(nline);
H.nx = str2double(nline{1});
H.ny = str2double(nline{2});
H.nz = str2double(nline{3});

% Define grid vectors. Tomography file must be defined on regular grid.
H.xg = H.xmin + H.dx*linspace(0,H.nx-1,H.nx);
H.yg = H.ymin + H.dy*linspace(0,H.ny-1,H.ny);
H.zg = H.zmin + H.dz*linspace(0,H.nz-1,H.nz);

% Fourth line of header: minimum/maximum values
char1 = '#';
while strcmp(char1,'#')
    nline = fgetl(fid);
    nline = strtrim(nline);
    char1 = nline(1);
    nhead = nhead + 1;
end
nline = strsplit(nline);
H.vpmin = str2double(nline{1});
H.vpmax = str2double(nline{2});
H.vsmin = str2double(nline{3});
H.vsmax = str2double(nline{4});
H.rmin = str2double(nline{5});
H.rmax = str2double(nline{6});

% Continuing reading until first data entry
char1 = '#';
while strcmp(char1,'#')
    nline = fgetl(fid);
    nline = strtrim(nline);
    if isempty(nline)
        char1 = '#';
    else
        char1 = nline(1);
    end
    nhead = nhead + 1;
end
nline = strsplit(nline);

% Check model file type
if length(nline) == 25
    tf_Cij = true;
elseif length(nline) == 6
    tf_Cij = false;
else
    error(['Expected 6 or 25 data columns but read ',num2str(length(nline)),'.']);
end
fclose(fid);

% Number of data lines
NMAX = H.nx*H.ny*H.nz;

% Read in model (dlmread is much faster than serial calls to fgetl)
U = dlmread(theFile,'',nhead,0);
U = U(:,4:end); % Delete coordinate columns (redunant information)

% Check number of data lines read
if ~(length(U(:,1)) == NMAX)
    error(['Expected to read ',num2str(NMAX),' lines but read ',num2str(length(U(:,1))),' lines.']);
end

% Define structure
if tf_Cij
    fldnms = {'RHO','C11','C12','C13','C14','C15','C16',...
                          'C22','C23','C24','C25','C26',...
                                'C33','C34','C35','C36',...
                                      'C44','C45','C46',...
                                            'C55','C56',...
                                                  'C66'};
else
    fldnms = {'Vp','Vs','RHO'};
end

% Define mode
for ii = 1:length(fldnms)
    data.(fldnms{ii}) = reshape(U(:,ii),H.nx,H.ny,H.nz);
end
clear U;

% Derive isotropic velocities
if tf_Cij
    data.Vp = (data.C11 + data.C22 + data.C33 + 2*data.C12 + 2*data.C13 + 2*data.C23)./9;
    data.Vs = (2*data.C11 + 2*data.C22 + 2*data.C33 + 6*data.C44 + 6*data.C55 + 6*data.C66 - 2*data.C12 - 2*data.C13 - 2*data.C23)./30;
    data.Vp = sqrt((10^9)*(data.Vp + (4/3)*data.Vs)./data.RHO);
    data.Vs = sqrt((10^9)*data.Vs./data.RHO);
    fldnms  = cat(2,fldnms,{'Vp','Vs'});
end

% Flip arrays?
if tf_flip
    for ii = 1:length(fldnms)
        data.(fldnms{ii}) = data.(fldnms{ii})(:,:,H.nz:-1:1);
    end
end

% Reset 'atmosphere' velocities?
if ~isempty(ireset)
    for n = 1:length(fldnms)
        for ii = 1:H.nx
            for jj = 1:H.ny
                ibad = find(squeeze(data.(fldnms{n})(ii,jj,:)) == ireset,1,'first');
                data.(fldnms{n})(ii,jj,ibad:end) = data.(fldnms{n})(ii,jj,ibad-1);
            end
        end
    end
end

% Flatten the model on Earth's curvature?
if tf_flat
    % Below interpolation introduces NaNs.
    % Manually extrapolate using center profile...
    iref = ceil(H.nx/2);
    jref = ceil(H.ny/2);
    
    R       = 6371000;
    [Y,X,Z] = meshgrid(H.yg,H.xg,H.zg);
    Z       = Z + sqrt((R^2) - (X.^2) - (Y.^2)) - R;
    for ii = 1:length(fldnms)
        data.(fldnms{ii}) = interp3(H.yg(:)',H.xg(:),H.zg(:)',data.(fldnms{ii}),Y,X,Z);
        % Extrapolate...
        ll = isnan(data.(fldnms{ii}));
        F  = squeeze(data.(fldnms{ii})(iref,jref,:));
        F  = permute(repmat(F(:),1,H.nx,H.ny),[2 3 1]);
        data.(fldnms{ii})(ll) = F(ll);
    end
end

% % Plot check results
% % Vp in XZ-plane
% figure;
% imagesc(H.xg./1000,H.zg./1000,squeeze(data.Vp(:,ceil(H.ny/2),:))');
% set(gca,'ydir','normal');
% pbaspect([(H.xmax-H.xmin)/(H.zmax-H.zmin), 1, 1]);
% colormap(jet);
% colorbar;
% title('XZ-Plane');
% xlabel('km');
% ylabel('km');
% 
% % Vp in YZ-plane
% figure;
% imagesc(H.yg./1000,H.zg./1000,squeeze(data.Vp(ceil(H.nx/2),:,:))');
% set(gca,'ydir','normal');
% pbaspect([(H.ymax-H.ymin)/(H.zmax-H.zmin), 1, 1]);
% colormap(jet);
% colorbar;
% title('YZ-Plane');
% xlabel('km');
% ylabel('km');
% 
% % Vp in XY-plane
% nz = 49;
% figure;
% imagesc(H.xg./1000,H.yg./1000,squeeze(data.Vp(:,:,nz))');
% set(gca,'ydir','normal');
% pbaspect([(H.xmax-H.xmin)/(H.ymax-H.ymin), 1, 1]);
% colormap(jet);
% colorbar;
% title('XY-Plane');
% xlabel('km');
% ylabel('km');
%
% % Magnitude of P-wave anisotropy in XY-plane at half-depth
% figure;
% imagesc(H.xg./1000,H.yg./1000,squeeze(Fp(:,:,ceil(H.nz/2)))');
% set(gca,'ydir','normal');
% pbaspect([(H.xmax-H.xmin)/(H.ymax-H.ymin), 1, 1]);
% colormap(jet);
% colorbar;
% title('XY-Plane');
% xlabel('km');
% ylabel('km');
% 
% % Orientation of P-wave anisotropy in XY-plane at half-depth
% % First, define anisotropy quivers
% n  = 2;
% s  = 250;
% izr = [31 H.nz];
% figure;
% for iz = izr(1):izr(2)
%     for ix = 1:n:H.nx
%         for iy = 1:n:H.ny
%             dx = s*abs(Fp(ix,iy,iz))*cos(phi(ix,iy,iz));
%             dy = s*abs(Fp(ix,iy,iz))*sin(phi(ix,iy,iz));
%             plot((H.xg(ix)/1000)+[dx,-dx],(H.yg(iy)/1000)+[dy,-dy],'-k'); hold on;
%         end
%     end
%     axis image;
%     box on;
%     title('XY-Plane Symmetry-Axis Orientation');
%     display(H.zg(iz)/1000);
%     shg;
%     pause;
%     clf;
% end
