%% Check Travel Time Table
% Simple script to compare tabled delays to analysis figures to make sure
% they are consistent.
close all;
clear;
clc;

% Input
evtid      = 10;
refModel   = 'ak135';
theEvent   = '~/research/CASCADIA/Waveforms_S/event_cascadia.mat';
theStation = '~/research/CASCADIA/Waveforms_S/station_cascadia.mat';
theATT     = '~/research/CASCADIA/Waveforms_S/ATT_cascadia.mat';
evtDir     = ['~/research/CASCADIA/Waveforms_S/PICKED/EVT',num2str(evtid)];
theCmap    = '';

% Load data structures
load(theEvent,'event');
load(theStation,'station');
load(theATT,'ATT');
% Plate boundaries
load('~/research/CASCADIA/P/m-files/PB_Cascadia.mat','PB');
if isempty(theCmap)
    Cmap = jet(64);
else
    load(theCmap,'Cmap');
end

%% Compute Predicted Travel Times
% Indexing
NARR     = length(ATT.arrivalTime);
[~,ievt] = ismember(evtid,event.id);
[~,ista] = ismember(ATT.station,station.name);
icalc    = ATT.event == evtid;

% Range
DLT = distance(station.latitude(ista),station.longitude(ista),...
    event.latitude(ievt),event.longitude(ievt));

% Compute 1D travel time predictions
TT1D = zeros(NARR,1);
for ii = 1:NARR
    if icalc(ii)
        TT1D(ii) = taup_time(refModel,ATT.phase{ii},DLT(ii),event.depth(ievt),0);
    end
end

% Measured travel times
TO = datenum(ATT.originTime);
TA = datenum(ATT.arrivalTime);
TT = (TA - TO)*24*60*60;

%% Prior Picks
CHN = {'T','T','PZ'};
WND = [0,10,10];

H = figure;
N = length(CHN);
DT_max = 0;
for ii = 1:N
    subplot(1,N,ii); hold on;
    plot(PB.boundaries(:,1),PB.boundaries(:,2),'-k','linewidth',2);
end
for ii = 1:N
    % Select Arrivals
    keep = (ATT.event == evtid) & strcmp(ATT.channel,CHN{ii}) & (ATT.window == WND(ii));
    
    % Unique arrivals check
    if ~(length(ATT.station(keep)) == length(unique(ATT.station(keep))))
        error('Non-unique selection!');
    end
    
    % Relative Delays
    DT = TT(keep) - TT1D(keep);
    DT = DT - mean(DT);
    % Track max absolute delay
    DT_max = max(max(abs(DT)),DT_max);
    
    % Plot delays
    subplot(1,N,ii);
    scatter(station.longitude(ista(keep)),station.latitude(ista(keep)),75,DT,'filled');
end
% Apply formatting
for ii = 1:N
    subplot(1,N,ii);
    box on; grid on; axis tight; axis image; colormap(Cmap);
    colorbar('southoutside'); caxis(DT_max*[-1,1]);
end
% Double figure width
H.Position(3) = 2*H.Position(3);

% Open analysis figure
eval(['!open ',evtDir,'/6_RelativeDelays.png']);

%% Plot Station-Averaged Delays

% All the original picks are flagged with 0 s window
iog = ATT.window == 0;
%iog = (ATT.window > 0) & strcmp(ATT.channel,'PZ');
N   = sum(iog);

% Reference travel-time structure
TTref.station = ATT.station(iog);
TTref.event   = ATT.event(iog);
TTref.phase   = ATT.phase(iog);
TTref.time    = zeros(N,1);

% Indexing
[~,ievt] = ismember(TTref.event,event.id);
[~,ista] = ismember(TTref.station,station.name);
% Compute distances
DLT = distance(event.latitude(ievt),event.longitude(ievt),station.latitude(ista),station.longitude(ista));
DPT = event.depth(ievt);

% Comptue TauP times 
for ii = 1:N
    TTref.time(ii) = taup_time(refModel,TTref.phase{ii},DLT(ii),DPT(ii),0);
    % If no phase found, what to do? How to interpret generic S phase?
    if isinf(TTref.time(ii)) && strcmp(TTref.phase{ii},'sS')
        TTref.time(ii) = taup_time(refModel,'sSdiff',DLT(ii),DPT(ii),0);
    end
end

% Original
% Measured travel times
TO = datenum(event.originDate(ievt));
TA = datenum(ATT.arrivalTime(iog));
TT = (TA - TO)*24*60*60;

% Initial Delays
DT_i = TT - TTref.time;

% Demeaning
MED = accumarray(ievt,DT_i,[],@mean,NaN);
RDT = DT_i - MED(ievt);
MSD = accumarray(ista,RDT,size(station.name),@mean,NaN);

% Plot
figure; hold on;
plot(PB.boundaries(:,1),PB.boundaries(:,2),'-k','linewidth',2);
scatter(station.longitude(ista),station.latitude(ista),50,MSD(ista),'filled');
axis image; colormap(Cmap); colorbar; caxis([-1.5,1.5]); box on; grid on;
