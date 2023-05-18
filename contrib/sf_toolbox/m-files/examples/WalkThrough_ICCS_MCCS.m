%% Walk Through Automatic Phase Picking
clear;
close all;
clc

% Required Paths
addpath(genpath('/Users/bvanderbeek/research/software/sf_toolbox/m-files'));
% Required Environment Vaiables
setenv('LOMAX','/Users/bvanderbeek/research/software/matlab_toolboxes/FilterPicker_Lomax');
setenv('SEGYMAT','/Users/bvanderbeek/research/software/matlab_toolboxes/segymat-master');

% Seismogram parameters
prjDir     = '/Users/bvanderbeek/research/NEWTON/AxiSEM_SPECFEM3D_Projects/projects/SubductionModel_101_Centered_Annulus_ANISO';
aPhase     = 'S';
evtid      = 1;
channel    = {'BXX','BXY'};
ext        = 'semv';
dataWin    = [-100 100]; % 300 to 1800 without moveout (start at 1200 for just S)
ttap       = 15;
filtopts   = [1/40 1/15 2 1]; % [f_min f_max order tf_zero_phase]
tf_rotate  = true;
tf_moveout = true;
aModel   = 'iasp91'; % 1D Model name

% Automatic picking parameters
ichan       = 1; % What channel index to pick
pick_window = [-50 50];
stw         = 200;
ltw         = 600;
ptw         = 50;
E1          = 10;
E2          = 10;
% ICCS parameters
iwstart = -5;
iwend   = 20;
itaper  = 15;
ccf_min = 0.7; % Reject any arrivals where the cross-correlation coefficient is less than this value
% MCC parameters
mwstart = -5;
mwend   = 15;
mtaper  = 15;

% Plotting parameters
n2plot   = 150;
sort_var = 'delta';
scale    = 5e6;
climits  = [-5 5];

%% Load data

% Stations
station = sf_read_stations(prjDir);
% % Remove stations?
% s2rmv = ((station.x < -800) | (station.x > 325)) | ((station.y < -1300) | (station.y > 1300));
% flds = fieldnames(station);
% for ii = 1:length(flds)
%     station.(flds{ii}) = station.(flds{ii})(~s2rmv);
% end

% Events
event = sf_read_events(prjDir);

% Seismograms
tic
data = sf_load_ascii_seis(station,event,evtid,channel,prjDir,ext,dataWin,...
    'filtopts',filtopts,'tf_rotate',tf_rotate,'tf_moveout',tf_moveout,'aPhase',aPhase,'aModel',aModel,'ttap',ttap);
toc

%data.baz = wrapTo180(data.baz);

% Index stations
[~,ist] = ismember(data.station,station.id);

% Make an initial pick of the seismograms--assume the 1D travel-time
% applied upon loading was sufficiently close
dtt = zeros(data.ntrace,1);

% Plot seismograms
sf_quick_plot_seis(ichan,data,scale,sort_var,dtt,n2plot);

%% Run ICCS alignment

% New ICCS routine
% % Lomax pick update
% clear PickParams;
% PickParams.method = 'lomax';
% PickParams.stw    = 30;
% PickParams.ltw    = 90;
% PickParams.ptw    = 10;
% PickParams.E1     = 10;
% PickParams.E2     = 10;
% PickParams.tadd   = 15;
% % STA/LTA pick update
% clear PickParams;
% PickParams.method = 'stalta';
% PickParams.stw    = 30;
% PickParams.ltw    = 30*3;
% PickParams.ttaper = 2;
% PickParams.pick_window = [-50 50];
% Manual Pick
clear PickParams
PickParams.method = 'manual';
PickParams.tadd = 0;
% % Constant pick
% clear PickParams;
% PickParams.method = 'nopick';
% PickParams.tpick  = tpick0;
[delays,tpick1,ibad,CCF,itol,S1] = ICCS(data.t(1),data.fs,squeeze(data.seis(:,ichan,:)),itaper,iwstart,iwend-iwstart,ccf_min,PickParams);
tpick1 = tpick1 + PickParams.tadd;
% Plot initial stacked trace and pick
figure;
plot(data.t,S1,'-b','linewidth',2);
hold on;
plot(tpick1*[1 1],rms(S1)*[-1 1],'-r','linewidth',2);
plot((tpick1+iwstart)*[1 1],rms(S1)*[-1 1],'--k','linewidth',1);
plot((tpick1+iwend)*[1 1],rms(S1)*[-1 1],'--k','linewidth',1);
box on;
grid on;
shg;

%% Look at correlation across array
figure;
scatter(station.longitude(ist),station.latitude(ist),40,CCF,'filled');
colormap(flipud(jet));
colorbar;
axis image;
box on;
title('Auto-correlation Normalized CCF');
caxis([0.5 1]);

% Dependence on back-azimuth and delta
figure;
plot3(data.baz,data.delta,CCF,'.b');
box on;
grid on;
xlabel('baz');
ylabel('delta');
title('CCF(baz,delta)');

%% Plot ICCS aligned seismograms and delays

H = sf_quick_plot_seis(ichan,data,scale,sort_var,delays+tpick1,n2plot);

% Plot MCC window
for ii = 1:length(H)
    figure(H{ii});
    plot(mwstart*[1 1],[H{ii}.Children.YTick(1),H{ii}.Children.YTick(end)],'-r','linewidth',2);
    plot(mwend*[1 1],[H{ii}.Children.YTick(1),H{ii}.Children.YTick(end)],'-r','linewidth',2);
end

figure;
scatter(station.longitude(ist),station.latitude(ist),40,delays,'filled');
colormap(jet);
colorbar;
axis image;
box on;
title('ICC Delay times');
caxis(climits);

%% Run MCC

[DT,S2] = MCCS(data.t(1),data.fs,squeeze(data.seis(:,ichan,:)),tpick1+delays,mtaper,mwstart,mwend-mwstart,CCF);

% Final Pick on MCC stack
tpick2 = tpick1;
% [tpick2,dtpick2] = filterPicker(1/data.fs,S2(:)',[stw ltw ptw E1 E2],1 + round(data.fs*(pick_window - data.t(1))));
% tpick2 = tpick2 + data.t(1);

figure;
%scatter(station.longitude(ist),station.latitude(ist),40,DT,'filled');
scatter(station.x(ist),station.y(ist),40,DT,'filled');
colormap(jet);
colorbar;
axis image;
box on;
title('MCC Delay times');
caxis(climits);
xlim([-1000 1000]);
ylim([-1500 1500]);


%% Save Picks?
% Save results
picks.eventid    = evtid;
picks.phase      = aPhase;
picks.channel    = data.channel(ichan);
picks.TRZ        = tf_rotate;
picks.station    = station.id;
picks.traveltime = tpick2 + data.tt1D + DT;
picks.filtopts   = filtopts;
% % Store lomax filter picking parameters
% picks.lomax.pick_window = pick_window;
% picks.lomax.stw         = stw;
% picks.lomax.ltw         = ltw;
% picks.lomax.ptw         = ptw;
% picks.lomax.E1          = E1;
% picks.lomax.E2          = E2;
% Store ICCS parameters
picks.iccs.iwstart  = iwstart;
picks.iccs.iwend    = iwend;
picks.iccs.itaper   = itaper;
picks.iccs.ccf_min  = ccf_min;
% Store MCCS parameters
picks.mccs.mwstart = mwstart;
picks.mccs.mwend   = mwend;
picks.mccs.mtaper  = mtaper;
save([prjDir,'/SEIS/CMTSOLUTION_',num2str(evtid),'/picks_',aPhase,'_',data.channel{ichan},'.mat'],'picks');

%% Final seismogram alignment

sf_quick_plot_seis(ichan,data,scale,sort_var,DT+tpick2,n2plot);

%% Look at frequency content of final stacked trace
% icenter = sqrt(station.x.^2 + station.y.^2);
% icenter = find(icenter == min(icenter));
% S       = squeeze(data.seis(strcmp(data.station,station.id(icenter)),ichan,:));
% S2      = S;

S = S2;
w = get_window(data.nsamp,data.fs,5,data.t(1),tpick2-5,25);
S = S(:).*w(:);

figure;
subplot(2,1,1); hold on;
plot(data.t,S,'-r','linewidth',2);
plot(data.t,S2,'-b');
xlim([-20 60])
box on; grid on;
xlabel('time (s)');
ylabel('amplitude');
title('Stacked Trace; Unfiltered');

fbin = linspace(0,data.fs,data.nsamp);
fbin = fbin - (fbin(end)/2);
Y2   = fftshift(fft(S2));
Y    = fftshift(fft(S));

subplot(2,1,2); hold on;
plot(fbin,abs(Y2),'-b');
plot(fbin,abs(Y),'-r','linewidth',2);
plot(ones(2,1)./15,[0,max(abs(Y))],'-k');
plot(ones(2,1)./13.1,[0,max(abs(Y))],'--k');
xlim([0 0.5]);
box on; grid on;
xlabel('frequency (Hz)');
ylabel('amplitude');
title('Frequency Content; Unfiltered');
