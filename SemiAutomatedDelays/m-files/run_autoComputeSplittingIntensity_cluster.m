% Run Automated Splitting Intensity
close all;
clear;
clc;

% Set environment variables
setenv('TAUPJAR',[getenv('STINGRAY'),'/contrib/TAUP/lib/TauP-2.4.5.jar']);
% Add required paths
addpath([getenv('STINGRAY'),'/toolbox/utils']);
addpath /storage2/unipd/vanderbeek/Cascadia_Joint/SplittingIntensity/m-files/mk_rd_mseed/

% Input
wrkDir     = '/storage2/unipd/vanderbeek/Cascadia_Joint/SplittingIntensity/data';
outFile    = [wrkDir,'/F_SplittingIntensity.mat']; % Output file with splitting intensity observations
dataDir    = [wrkDir,'/waveforms']; % Waveform directory
theATT     = [wrkDir,'/ATT_cascadia.mat']; % The arrival time table
theEvent   = [wrkDir,'/event_cascadia.mat']; % Event structure
theStation = [wrkDir,'/station_cascadia.mat']; % Station structure
refModel   = 'ak135'; % Reference 1D model

% Load data
load(theATT,'ATT');
load(theEvent,'event');
load(theStation,'station');

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
npool = 6;
% Start parallel pool
parpool('local',npool,'AttachedFiles',{getenv('TAUPJAR')});
tic
parfor (ii = 1:nevt, npool)
    evtid = strsplit(EVTX(ii).name,'EVT');
    evtid = str2double(evtid{end});
    F(ii) = autoComputeSplittingIntensity(evtid,ATT,TT1D,event,station,dataDir,refModel);
    fprintf(['\n Finished event ',num2str(ii),' of ',num2str(nevt),'. \n']);
end
toc
% Delete parallel pool
delete(gcp);
save(outFile,'F');
