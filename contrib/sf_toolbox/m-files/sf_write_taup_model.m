function [z,vp,vs,rho] = sf_write_taup_model(prjDir)

% Read AxiSEM model
[R,rho,Vph,Vpv,Vsh,Vsv] = sf_read_dsm_model([prjDir,'/MESH/iasp91_dsm']);

% Depth
z = 6371 - R;
z = flipud([z(:,2),z(:,1)]);
z = mean(z,2);
z = [0; z; 6371];

% Vp (km/s)
vp = flipud((Vph(:,1) + Vpv(:,1))./2);
vp = [vp(1); vp; vp(end)];
% Vs (km/s)
vs = flipud((Vsh(:,1) + Vsv(:,1))./2);
vs = [vs(1); vs; vs(end)];
% Density (g/m^3)
rho = flipud(rho);
rho = [rho(1); rho; rho(end)];

% Identify boundaries
ioc = find(vs==0,1,'first'); % Outer core
iic = find(vs==0,1,'last'); % Inner core

% Check/create file
if isfile([prjDir,'/model1D.nd'])
    error('The file model1D.nd already exists!');
else
    fid = fopen([prjDir,'/model1D.nd'],'w');
end

for ii = 1:length(vp)
    if ii == ioc
        fprintf(fid,'%10.3f %8.3f %8.3f %8.3f\n',z(ii),vp(ii-1),vs(ii-1),rho(ii-1));
        fprintf(fid,'%s\n','outer-core');
        fprintf(fid,'%10.3f %8.3f %8.3f %8.3f\n',z(ii),vp(ii),vs(ii),rho(ii));
    end
    fprintf(fid,'%10.3f %8.3f %8.3f %8.3f\n',z(ii),vp(ii),vs(ii),rho(ii));
    if ii == iic
        fprintf(fid,'%10.3f %8.3f %8.3f %8.3f\n',z(ii),vp(ii),vs(ii),rho(ii));
        fprintf(fid,'%s\n','inner-core');
        fprintf(fid,'%10.3f %8.3f %8.3f %8.3f\n',z(ii),vp(ii+1),vs(ii+1),rho(ii+1));
    end
end


% This created discontinuities at every interface
% % Depth (km)
% z   = mean(6371 - R,2);
% z   = flipud([6371; z; 0]);
% % Vp
% vp  = (Vph(:,1) + Vpv(:,1))./2;
% vp  = flipud([vp(1); vp; vp(end)]);
% % Vs
% vs  = (Vsh(:,1) + Vsv(:,1))./2;
% vs  = flipud([vs(1); vs; vs(end)]);
% % Density (g/cm
% rho = rho(:,1);
% rho = flipud([rho(1); rho; rho(end)]);
% 
% % Check/create file
% if isfile([prjDir,'/model1D.taup'])
%     error(['The file ',theFile,' already exists!']);
% else
%     fid = fopen([prjDir,'/model1D.tvel'],'w');
% end
% % Two header lines
% fprintf(fid,'%s\n','#');
% fprintf(fid,'%s\n','#');
% for ii = 1:length(z)
%     fprintf(fid,'%10.3f %8.3f %8.3f %8.3f\n',z(ii),vp(ii),vs(ii),rho(ii));
% end
% fclose(fid);
