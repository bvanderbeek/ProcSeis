%% Input
% Must define path to TauP jar file
setenv('TAUPJAR','../lib/TauP-2.4.5.jar');

aModel = 'iasp91'; % A TauP-recognized 1D model name (see TauP manual)
aPhase = 'P'; % A TauP-recognized seismic phase (see TauP manual)
delta  = 80; % Source-receiver distance in arc degrees
sdepth = 50; % Source depth; positive into the Earth (km)
rdepth = 0; % Receiver depth; positive into the Earth (km)

%% Compute a Traveltime
[tt,rayP,inc,toff,outPhase,iwarn] = taup_time(aModel,aPhase,delta,sdepth,rdepth);

%% Compute a Ray Path
% Define geographic source-receiver positions
slat = 0;
slon = 0;
rlat = 0;
rlon = delta;
% Get ray path
[latray,lonray,zray,ttray,dray,outPhase] = taup_path(aModel,aPhase,slat,slon,sdepth,rlat,rlon);
% Plot ray path
Re   = 6371; % Earth radius
rray = Re + zray;
xray = rray.*sind(dray);
yray = rray.*cosd(dray);
figure(101); hold on;
plot(Re*cosd(linspace(0,360,361)),Re*sind(linspace(0,360,361)),'-k','linewidth',2);
plot(xray,yray,'-r','linewidth',2);
plot(xray(1),yray(1),'*r','markersize',10);
plot(xray(end),yray(end),'vb','markersize',10);
axis image; box on; grid on;

%% Load a 1D Radial Velocity Model
% Define depth vector; positive into the Earth (km)
Rz = linspace(0,1000,1001)';
% Get 1D model
[vp,vs,rho,Qp,Qs,Rz] = get_TauP_model(Rz,aModel);
% Plot velocity
figure(201); hold on;
plot(vp,-Rz,'-b','linewidth',2);
plot(vs,-Rz,'-r','linewidth',2);
pbaspect([1,2,1]); box on; grid on;
xlabel('velocity (km/s)');
ylabel('depth (km)');
legend('Vp','Vs');
