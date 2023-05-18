% Run Automated Splitting Intensity
close all;
clear;
clc;

% Set environment variables
setenv('TAUPJAR','/Users/bvanderbeek/research/software/TauP-2.4.5/lib/TauP-2.4.5.jar');
% Add required paths
addpath ~/research/software/Stingray_GIT_BPV/toolbox/utils/
addpath ../../codes/mk_rd_mseed/

% Input
outFile    = 'F_SplittingIntensity.mat'; % Output file with splitting intensity observations
dataDir    = '~/research/CASCADIA/Waveforms_S/PICKED'; % Waveform directory
theATT     = '~/research/CASCADIA/Waveforms_S/ATT_cascadia.mat'; % The arrival time table
theEvent   = '~/research/CASCADIA/Waveforms_S/event_cascadia.mat'; % Event structure
theStation = '~/research/CASCADIA/Waveforms_S/station_cascadia.mat'; % Station structure
refModel   = 'ak135'; % Reference 1D model

% Load data
load(theATT,'ATT');
load(theEvent,'event');
load(theStation,'station');

% Load plate boundaries
load('~/research/CASCADIA/Joint/m-files/PB_Cascadia.mat','PB');

%% Compute Reference Travel-Times
% Source-receiver distance for every arrival
[~,ista] = ismember(ATT.station,station.name);
[~,ievt] = ismember(ATT.event,event.id);
DLT      = distance(station.latitude(ista),station.longitude(ista),event.latitude(ievt),event.longitude(ievt));

% Predicted travel-times
nobs = length(ATT.station);
TT1D = zeros(nobs,1);
for ii = 1:nobs
    TT1D(ii) = taup_time(refModel,ATT.phase{ii},DLT(ii),event.depth(ievt(ii)),0);
end

%% Compute Splitting Intensity
% Get list of events
EVTX = dir([dataDir,'/EVT*']);
nevt = length(EVTX);
% Initialize splitting intensity structure
F(nevt).eventid = [];
F(nevt).station = [];
F(nevt).dlt     = [];
F(nevt).baz     = [];
F(nevt).paz     = [];
F(nevt).SI      = [];
F(nevt).SIj     = [];
% Loop over events
npool = 0;
tic
parfor (ii = 1:nevt, npool)
    evtid = strsplit(EVTX(ii).name,'EVT');
    evtid = str2double(evtid{end});
    F(ii) = autoComputeSplittingIntensity(evtid,ATT,TT1D,event,station,dataDir,refModel);
    fprintf(['\n Finished event ',num2str(ii),' of ',num2str(nevt),'. \n']);
end
toc
save(outFile,'F');

%% Fit Splitting Intensity
% Load splitting intensity structure
% *_inorm1, *_inorm1_pz, *_inorm2_pz, *_inorm3_pz
load('F_SplittingIntensity_inorm1.mat','F');
% Demean splitting intensity observations
tf_demean = true;
% Concatonate structure fields
G.eventid = cat(1,F(:).eventid);
G.station = cat(1,F(:).station);
G.dlt = cat(1,F(:).dlt);
G.baz = cat(1,F(:).baz);
G.paz = cat(1,F(:).paz);
G.SI = cat(1,F(:).SI);
G.SIj = cat(1,F(:).SIj);

% Define data to fit
SI  = G.SIj(:,1);
dSI = G.SIj(:,2);
% Option to demean by event
if tf_demean
    EID = unique(G.eventid);
    [~,ievt] = ismember(G.eventid,EID);
    MEO = accumarray(ievt,SI,[],@mean);
    SI = SI - MEO(ievt);
end

% Define observation azimuth positive counter-clockwise from east
% Backazimuth is positive clockwise of north
RAZ = 270 - G.baz;
% Define polarization angle positive counter-clockwise from east
% Stored polarisation azimuth is postive counter-clockwise from transverse
PAZ = G.paz - 90;
% Observation azimuth is
OAZ = RAZ + PAZ;

% Fit station splitting intensity
keep = (abs(SI) < 3) & (abs(SI./dSI) > 1.5);
F = fit_splitting_intensity(SI(keep),OAZ(keep),G.station(keep),dSI(keep));

% Station Indexing
[~,ista] = ismember(F.UID,station.name);

% NaN bad fits
ibad = abs(F.DT_2) > 4;
F.DT_2(ibad) = NaN;
F.DT_1(ibad) = NaN;
F.DT_0(ibad) = NaN;

% Plot results
s    = 0.2;
cint = 0:0.25:3;
plot_splits(station.longitude(ista),station.latitude(ista),F.DT_2,F.AZM_2,cint,s);
% Plot plate boundaries
plot(PB.boundaries(:,1),PB.boundaries(:,2),'-k','linewidth',2);
title(num2str(sum(keep)));

%% Fit Geographic Binned Data
% Grid spacing adjusted for latitude
dlat = 0.1;
dlon = distance(mean(station.latitude),mean(station.longitude) - 0.5,...
    mean(station.latitude),mean(station.longitude) + 0.5);
dlon = dlat/dlon;

% Grid extents
minlon = floor(min(station.longitude(:)));
maxlon = ceil(max(station.longitude(:)));
minlat = floor(min(station.latitude(:)));
maxlat = ceil(max(station.latitude(:)));
% Number of points
Nlon = 1 + round((maxlon - minlon)/dlon);
Nlat = 1 + round((maxlat - minlat)/dlat);
% Grid vectors
glon = linspace(minlon,maxlon,Nlon);
glat = linspace(minlat,maxlat,Nlat);
ind  = reshape(1:(Nlon*Nlat),Nlon,Nlat);

% Get coordinates for each split
[~,ista] = ismember(G.station,station.name);
LAT = station.latitude(ista);
LON = station.longitude(ista);

% Associate each split to a grid point
igrd = interp2(glat(:)',glon(:),ind,LAT,LON,'nearest');
igrd = round(igrd);
% Construct linear index grid
[glat,glon] = meshgrid(glat,glon);
glat = glat(:);
glon = glon(:);

% Fit station splitting intensity
F = fit_splitting_intensity(SI,OAZ,igrd,dSI);

% NaN bad fits
ibad = abs(F.DT_2) > 4;
F.DT_2(ibad) = NaN;
F.DT_1(ibad) = NaN;
F.DT_0(ibad) = NaN;

% Plot results
s    = 0.2;
cint = 0:0.25:3;
plot_splits(glon(F.UID),glat(F.UID),F.DT_2,F.AZM_2,cint,s);
% Plot plate boundaries
plot(PB.boundaries(:,1),PB.boundaries(:,2),'-k','linewidth',2);

%% Observed Cascadia Splits
load('S_WUS_Splits.mat','S');

% Plot results
s    = 0.2;
cint = 0:0.25:3;
H = plot_splits(S.longitude,S.latitude,S.delay,90-S.azimuth,cint,s);
% Plot plate boundaries
plot(PB.boundaries(:,1),PB.boundaries(:,2),'-k','linewidth',2);

%% Ad-hoc add SI measurements
% Change plotting parameters in function
%plot_splits(glon(F.UID),glat(F.UID),F.DT_2,F.AZM_2,cint,s,H);

plot_splits(station.longitude(ista),station.latitude(ista),F.DT_2,F.AZM_2,cint,s,H);
%% Function: Fit Splitting Intensity
function F = fit_splitting_intensity(SI,OAZ,ID,varargin)
if length(varargin) == 1
    w = varargin{1};
else
    w = ones(size(SI));
end
% Unique location IDs corresponding to splitting intensity observations
F.UID = unique(ID);
NID = length(F.UID);
% Allocate storage arrays
F.DT_2  = NaN(NID,1);
F.AZM_2 = NaN(NID,1);
F.DT_1  = NaN(NID,1);
F.AZM_1 = NaN(NID,1);
F.DT_0  = NaN(NID,1);
F.NOBS  = NaN(NID,1);
for ii = 1:NID
    % Select location
    k = ismember(ID,F.UID(ii));
    F.NOBS(ii) = sum(k);
    % Fit sinusoidal curves
    % sin(a - b) = sin(a)*cos(b) - cos(a)*sin(b)
    % a = psi = symmetry axis azimuth
    % b = oaz = observations azimuth
    % Hence, the negative in the sine-terms below
    if F.NOBS(ii) > 10
        % Solve system
        A = [ones(F.NOBS(ii),1),-sind(2*OAZ(k)),cosd(2*OAZ(k)),-sind(OAZ(k)),cosd(OAZ(k))];
        %x = lsqr(A,SI(k),[],1000);
        x = lscov(A,SI(k),w(k));
        % Splitting function is sin(2theta) NOT cos(2theta). Hence, the
        % sign change and 90 phase shift with respect to what we normally
        % use for cos(2theta) parameterization.
        F.AZM_2(ii) = atan2d(x(3),x(2))/2;
        F.DT_2(ii)  = sqrt(sum(x(2:3).^2));
        % The 1-theta terms
        F.AZM_1(ii) = atan2d(x(5),x(4));
        F.DT_1(ii)  = sqrt(sum(x(4:5).^2));
        % The mean offset
        F.DT_0(ii) = x(1);
    elseif (F.NOBS(ii) > 4) && (F.NOBS(ii) < 10)
        % Solve system
        A = [ones(F.NOBS(ii),1),-sind(2*OAZ(k)),cosd(2*OAZ(k)),-sind(OAZ(k)),cosd(OAZ(k))];
        %x = lsqr(A,SI(k),[],1000);
        x = lscov(A,SI(k),w(k));
        % Splitting function is sin(2theta) NOT cos(2theta). Hence, the
        % sign change and 90 phase shift with respect to what we normally
        % use for cos(2theta) parameterization.
        F.AZM_2(ii) = atan2d(x(3),x(2))/2;
        F.DT_2(ii)  = sqrt(sum(x(2:3).^2));
        % The mean offset
        F.DT_0(ii) = x(1);
    end
%     if F.NOBS(ii) > 50
%         xp  = -180:180;
%         SIp = x(1) + x(2)*sind(2*xp) - x(3)*cosd(2*xp) + x(4)*sind(xp) - x(5)*cosd(xp);
%         H = figure; hold on;
%         plot(wrapTo180(OAZ(k)),SI(k),'.k','markersize',25);
%         plot(xp,SIp,'-r','linewidth',2);
%         box on; grid on; axis tight; axis square; shg;
%         title(num2str(F.NOBS(ii)));
%         pause;
%         close(H);
%     end
end
end

%% Function: Plot Splits
function H = plot_splits(x,y,DELAY,AZIM,cint,s,varargin)
if length(varargin) == 1
    H = varargin{1};
else
    H = figure;
end
% Create scaled quivers
dx = s*DELAY.*cosd(AZIM);
dy = s*DELAY.*sind(AZIM);
dx = x + [-dx,dx];
dy = y + [-dy,dy];

% Plot color-coded qivers
figure(H); hold on;
dC   = cint(end)/(length(cint)-1);
Cmap = colormap(parula(length(cint)-1));
for ista = 1:length(DELAY)
    icolor = floor(DELAY(ista)/dC);
    icolor = max(icolor,1);
    icolor = min(icolor,length(cint)-1);
    %plot(dx(ista,:)',dy(ista,:)','-','linewidth',2,'Color',Cmap(icolor,:));
    plot(dx(ista,:)',dy(ista,:)','-r','linewidth',1);
end
axis image; box on; grid on; axis tight;

% Construct colorbar
cb = colorbar('southoutside');
cb.Ticks = linspace(0,1,length(cint));
cb.TickLabels = cellstr(num2str(cint(:)));

end