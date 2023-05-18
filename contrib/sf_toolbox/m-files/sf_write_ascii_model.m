function sf_write_ascii_model(outFile,nx,xmin,xmax,ny,ymin,ymax,nz,zmin,zmax,tf_Cij,varargin)
% Write out ASCII model to be read by SPECFEM3D
%
% Automatically sets min and max values to 95% and 105% of those already
% present
%
% This routine and SPECFEM3D assume Cij is in GPa and density is in kg/m^3!!!
%
% Beware of SPECFEM3D maximum string length (512). Current format produces
% data rows that are ~300 characters.

% Do not let Vp, Vs, or RHO extend beyond these values
vpmin  = 0; vpmax  = 50000;
vsmin  = 0; vsmax  = 50000;
rhomin = 0; rhomax = 50000;

if length(varargin) == 3
    % Isotropic
    RHO = varargin{1};
    a11 = varargin{2};
    b44 = varargin{3};
    a33 = a11;
    b66 = b44;
    eta = 1;
    dip = []; % Used a check for rotation; no rotation if empty
elseif length(varargin) == 8
    % Hexagonal symmetry
    RHO = varargin{1};
    a11 = varargin{2};
    a33 = varargin{3};
    b44 = varargin{4};
    b66 = varargin{5};
    eta = varargin{6};
    azim = varargin{7};
    dip  = varargin{8};
elseif length(varargin) == 22
    % Custom
    RHO  = varargin{1};
    C11  = varargin{2};
    C12  = varargin{3};
    C13  = varargin{4};
    C14  = varargin{5};
    C15  = varargin{6};
    C16  = varargin{7};
    C22  = varargin{8};
    C23  = varargin{9};
    C24  = varargin{10};
    C25  = varargin{11};
    C26  = varargin{12};
    C33  = varargin{13};
    C34  = varargin{14};
    C35  = varargin{15};
    C36  = varargin{16};
    C44  = varargin{17};
    C45  = varargin{18};
    C46  = varargin{19};
    C55  = varargin{20};
    C56  = varargin{21};
    C66  = varargin{22};
else
    error('Incorrect number of variable input arguments')
end

% Define 21 elastic constants if needed
if ((tf_Cij) && (length(varargin) == 3)) || (length(varargin) == 8)
    % Elastic constants for isotropic or hexagonal symmetry
    C11  = RHO.*a11.*a11./(10^9); % Symmetry plane P-wave velocity
    C22  = RHO.*a11.*a11./(10^9); % Symmetry plane P-wave velocity
    C33  = RHO.*a33.*a33./(10^9); % Symmetry AXIS P-wave velocity
    C44  = RHO.*b44.*b44./(10^9); % Symmetry plane S-wave velocity
    C55  = RHO.*b44.*b44./(10^9); % Symmetry plane S-wave velocity
    C66  = RHO.*b66.*b66./(10^9); % Symmetry AXIS P-wave velocity
    C12  = (C11 - 2*C66);
    C13  = eta.*(C11 - 2*C44);
    C23  = C13;
    % All other components are zero valued
    C14  = zeros(size(RHO));
    C15  = zeros(size(RHO));
    C16  = zeros(size(RHO));
    C24  = zeros(size(RHO));
    C25  = zeros(size(RHO));
    C26  = zeros(size(RHO));
    C34  = zeros(size(RHO));
    C35  = zeros(size(RHO));
    C36  = zeros(size(RHO));
    C45  = zeros(size(RHO));
    C46  = zeros(size(RHO));
    C56  = zeros(size(RHO));
    
    % Rotate tensors
    if (length(varargin) == 8) && ~isempty(dip)
        if (length(dip) == 1) && isempty(azim)
            disp('Rotating symmetry axis...');
            % Diagonals
            Csym = C33;
            C33  = C11;
            C11  = Csym;
            Csym = C66;
            C66  = C44;
            C44  = Csym;
            % Off-diagonals
            Csym = C12;
            C12  = C23;
            C23  = Csym;
            clear Csym
        else
            disp('Rotating hexagonal Cij tensors...');
            for n = 1:(nx*ny*nz)
                Cij = [C11(n) C12(n) C13(n) C14(n) C15(n) C16(n);...
                    C12(n) C22(n) C23(n) C24(n) C25(n) C26(n);...
                    C13(n) C23(n) C33(n) C34(n) C35(n) C36(n);...
                    C14(n) C24(n) C34(n) C44(n) C45(n) C46(n);...
                    C15(n) C25(n) C35(n) C45(n) C55(n) C56(n);...
                    C16(n) C26(n) C36(n) C46(n) C56(n) C66(n)];
                Cij = MS_rot3(Cij,0,dip(n),azim(n)); %,'order',[3 2 1]);
                % Store rotated values
                C11(n) = Cij(1,1);
                C12(n) = Cij(1,2);
                C13(n) = Cij(1,3);
                C14(n) = Cij(1,4);
                C15(n) = Cij(1,5);
                C16(n) = Cij(1,6);
                C22(n) = Cij(2,2);
                C23(n) = Cij(2,3);
                C24(n) = Cij(2,4);
                C25(n) = Cij(2,5);
                C26(n) = Cij(2,6);
                C33(n) = Cij(3,3);
                C34(n) = Cij(3,4);
                C35(n) = Cij(3,5);
                C36(n) = Cij(3,6);
                C44(n) = Cij(4,4);
                C45(n) = Cij(4,5);
                C46(n) = Cij(4,6);
                C55(n) = Cij(5,5);
                C56(n) = Cij(5,6);
                C66(n) = Cij(6,6);
            end
        end
    end
end

% Derive average velocities
if (length(varargin) == 8) || (length(varargin) == 22)
    Vp = (C11 + C22 + C33 + 2*C12 + 2*C13 + 2*C23)./9;
    Vs = (2*C11 + 2*C22 + 2*C33 + 6*C44 + 6*C55 + 6*C66 - 2*C12 - 2*C13 - 2*C23)./30;
    Vp = sqrt((10^9)*(Vp + (4/3)*Vs)./RHO);
    Vs = sqrt((10^9)*Vs./RHO);
else
    Vp = a11;
    Vs = b44;
end
clear a11 a33 b44 b66 mkiso

% Grid spacing
dx = (xmax-xmin)/(nx-1);
dy = (ymax-ymin)/(ny-1);
dz = (zmax-zmin)/(nz-1);

% Coordinate arrays
X = xmin + dx*linspace(0,nx-1,nx);
Y = ymin + dy*linspace(0,ny-1,ny);
Z = zmin + dz*linspace(0,nz-1,nz); % Assumes first plane of nodes is at depth
[Y,X,Z] = meshgrid(Y,X,Z); % Assumes x-coordinate follows first dimension (i.e. rows)
% The wrong way to define coordinates
% [X,Y,Z] = meshgrid(X,Y,Z); % Assumes x-coordinate follows second dimension (i.e. columns)
% Z = zmax - dz*linspace(0,nz-1,nz); % Assumes first plane of nodes is at surface

% Open file
% I think to use the internal mesher (xmeshfem3D) the model file must be
% named 'tomography_model.xyz'
if exist(outFile,'file') > 0
    error(['The file ', outFile, 'already exists.']);
else
    % I think to use the internal mesher (xmeshfem3D) the model file must
    % be named 'tomography_model.xyz'
    fid = fopen(outFile,'w');
end

% Write model boundary header (X_MIN Y_MIN Z_MIN X_MAX Y_MAX Z_MAX)
fprintf(fid,'%1s %10s %12s %12s %12s %12s %12s\n',...
    '#','X_MIN','Y_MIN','Z_MIN','X_MAX','Y_MAX','Z_MAX');
fprintf(fid,'%12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n',...
    xmin,ymin,zmin,xmax,ymax,zmax);

% Write model spacing header (DX DY DZ)
fprintf(fid,'%1s %10s %12s %12s\n','#','DX','DY','DZ');
fprintf(fid,'%12.3f %12.3f %12.3f\n',dx,dy,dz);

% Write model dimensions header (NX NY NZ)
fprintf(fid,'%1s %5s %7s %7s\n','#','NX','NY','NZ');
fprintf(fid,'%7.0f %7.0f %7.0f\n',nx,ny,nz);

% Write model value range header (VP_MIN VP_MAX VS_MIN VS_MAX RHO_MIN RHO_MAX...QP? QS?)
ivp =  (Vp > vpmin) & (Vp < vpmax); if ~any(ivp(:)); ivp(1) = true; end
ivs  = (Vs > vsmin) & (Vs < vsmax); if ~any(ivs(:)); ivs(1) = true; end
irho = (RHO > rhomin) & (RHO < rhomax); if ~any(irho(:)); irho(1) = true; end
fprintf(fid,'%1s %6s %8s %8s %8s %8s %8s\n',...
    '#','VP_MIN','VP_MAX','VS_MIN','VS_MAX','RHO_MIN','RHO_MAX');
fprintf(fid,'%8.1f %8.1f %8.1f %8.1f %8.1f %8.1f\n',...
    floor(0.95*min(Vp(ivp))),ceil(1.05*max(Vp(ivp))),floor(0.95*min(Vs(ivs))),ceil(1.05*max(Vs(ivs))),...
    floor(0.95*min(RHO(irho))),ceil(1.05*max(RHO(irho))));

% Write the model
if tf_Cij
    fprintf(fid,['\n%1s %10s %12s %12s %8s %12s %12s %12s %12s %12s %12s %12s ',...
                                          '%12s %12s %12s %12s %12s %12s %12s ',...
                                          '%12s %12s %12s %12s %12s %12s %12s\n'],...
        '#','X','Y','Z','RHO','C11','C12','C13','C14','C15','C16',...
                                    'C22','C23','C24','C25','C26',...
                                          'C33','C34','C35','C36',...
                                                'C44','C45','C46',...
                                                      'C55','C56',...
                                                            'C66');
    % Write model
    fprintf(fid,['%12.3f %12.3f %12.3f %8.1f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f ',...
                                            '%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f ',...
                                            '%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n'],...
        [X(:),Y(:),Z(:),RHO(:),C11(:),C12(:),C13(:),C14(:),C15(:),C16(:),...
                                      C22(:),C23(:),C24(:),C25(:),C26(:),...
                                              C33(:),C34(:),C35(:),C36(:),...
                                                     C44(:),C45(:),C46(:),...
                                                            C55(:),C56(:),...
                                                                   C66(:)]');
    disp('Wrote a Cij model file.');
else
    fprintf(fid,'\n%1s %10s %12s %12s %8s %8s %8s\n',...
        '#','X','Y','Z','VP','VS','RHO');
    fprintf(fid,'%12.3f %12.3f %12.3f %8.1f %8.1f %8.1f\n',...
        [X(:),Y(:),Z(:),Vp(:),Vs(:),RHO(:)]');
    disp('Wrote an isotropic model file.');
end

% Close file
fclose(fid);
