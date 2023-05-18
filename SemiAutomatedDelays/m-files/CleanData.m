%% Data Cleaning
% Re-organize station, event, and arrival time data into cleaner
% structures with my own conventions for ease of use. Should only have to
% run this once to make clean structures from Miles' data.
close all;
clear;
clc;

% The master data file from Miles with the clPick, srStation, and srEvent
% structures.
theData = '../data/cascadia_m58-All_years_0-90_S.mat';
% The horizontal channel orientation information
theHorizontals = '../data/orientations/horizontals_master_DEC2022.mat';

% Load data
load(theData,'srEvent','srStation','clPick');

% Load horizontal orientations
load(theHorizontals,'horizontals');

%% Create New Station File
station.network   = srStation.network(:);
station.locid     = srStation.locationcode(:);
station.name      = srStation.name(:);
station.longitude = srStation.longitude(:);
station.latitude  = srStation.latitude(:);
station.elevation = srStation.elevation(:);
station.startDate = srStation.start_date(:);
station.endDate   = srStation.end_date(:);
% Store channel orientations
[lsta,ista] = ismember(station.name,horizontals.station);
if (sum(lsta) == length(lsta))
    station.channel         = horizontals.channel;
    station.channel_azimuth = horizontals.channel_azimuth(ista,:);
else
    error('Missing channel information!');
end
% Clock Drift
station.clockDrift = srStation.clock_drift(:);
% Number of stations
station.nsta = length(station.name);

% Convert dates to Matlab's datetime format. Note that the 'format'
% specifier only refers to the display format and not the actual precision
% of the date data.
station.startDate = datetime(station.startDate,'InputFormat','yyyy-MM-dd HH:mm:ss','format','yyyy-MM-dd HH:mm:ss.SSS');
station.endDate   = datetime(station.endDate,'InputFormat','yyyy-MM-dd HH:mm:ss','format','yyyy-MM-dd HH:mm:ss.SSS');

% Check for unique names
if ~(length(unique(station.name)) == station.nsta)
    error('Non-unique station names!');
end

%% Create New Event Structure
event.id          = srEvent.id(:);
event.longitude   = srEvent.longitude(:);
event.latitude    = srEvent.latitude(:);
event.depth       = srEvent.depth(:);
event.originDate  = datetime(srEvent.origintime(:),'ConvertFrom','posixtime','format','yyyy-MM-dd HH:mm:ss.SSS');
event.Mw          = srEvent.magnitude(:);
event.nevt        = length(event.id);

% Check for unique events
if ~(length(unique(event.id)) == event.nevt)
    error('Non-unique event identifiers!');
end

%% Check for Unique Arrivals
% Check for non-unique data
[~,ista]  = ismember(clPick.station,station.name);
[~,ievt]  = ismember(clPick.eventid,event.id);
[~,iphs]  = ismember(clPick.phase,unique(clPick.phase));
[~,iu,id] = unique([ievt(:),ista(:),iphs(:)],'rows','stable');
% Check duplicate differences
NDUP = accumarray(id,ones(size(id)),[],@sum)'; % Number of duplicates
DAT  = accumarray(id,clPick.measured_arrival_time,[],@range)'; % Maximum difference

fprintf(['\n ',num2str(sum(NDUP > 1)),' arrivals have duplicate entries. \n']);

% None of the arrival time differences seem to be too extreme so let's just
% average duplicate arrival times
AT  = accumarray(id,clPick.measured_arrival_time,[],@mean)'; % Mean arrival
MRE = accumarray(id,clPick.residual_error,[],@mean)'; % Mean residual error

% Remove duplicate entries
flds = fieldnames(clPick);
for ii = 1:length(flds)
    clPick.(flds{ii}) = clPick.(flds{ii})(iu);
end
% Define averaged values
clPick.measured_arrival_time = AT;
clPick.residual_error        = MRE;
% Redefine arrival time residual
clPick.residual = clPick.measured_arrival_time - clPick.predicted_arrival_time;

%% Build Arrival Time Table
% Variable Descriptions
%       Event: Unique numeric event identifier
%  originTime: Event origin time (Matlab datetime format)
%     station: Unique character string identifying station
%     channel: Unique character string describing channel
%       phase: Unique character string descringing seismic phase picked
%         pol: Polarisation [azimuth,elevation] measured in ray-aligned
%              TQL-coordinates (deg.; positive-counter-clockwise).
%      filter: A 1x5 vector describing filter parameters
%              1) Filter type; 0 for Butterworth
%              2) Minimum corner frequency (Hz)
%              3) Maximum corner frequency (Hz)
%              4) Order of filter
%              5) Zero phase filter (0) or minimum phase (1)
%      window: Length of cross-correlation window (s). Window starts at
%               measured arrival time.
% arrivalTime: Measured arrival time (Matlab datetime format)
%       error: Measurement error (s)
%  lastUpdate: Last time measurement was modified (Matlab datetime format)

% Table Name
theTable = '../data/ATT_cascadia.mat';

% Table column names and types
varNames = {'event','originTime','network','station','channel','phase',...
    'pol','filter','window','arrivalTime','error','lastUpdate'};
varTypes = {'double','datetime','cellstr','cellstr','cellstr','cellstr',...
    'double','double','double','datetime','double','datetime'};

% Initialize empty table
% ATT = table('VariableNames',varNames,'Size',size(varNames),'VariableTypes',varTypes);

% Build table from clPick structure
[~,ista] = ismember(clPick.station(:),station.name);
[~,ievt] = ismember(clPick.eventid(:),event.id);
nar = length(clPick.measured_arrival_time);
OT  = event.originDate(ievt);
NTW = station.network(ista);
CHN = repmat({'T'},nar,1); % ASSUMED CHANNEL
POL = zeros(nar,2); % ASSUMED POLARISATION AZIMUTH,ELEVATION
FLT = repmat([0,1/33,1/12,3,0],nar,1); % ASSUMED FILTER
WND = repmat(15,nar,1); % ASSUMED WINDOW
AT  = datetime(clPick.measured_arrival_time(:),'ConvertFrom','posixtime','format','yyyy-MM-dd HH:mm:ss.SSS');
LUD = repmat(datetime,nar,1);
% Make table
ATT = table(clPick.eventid(:),OT,NTW,clPick.station(:),CHN,clPick.phase(:),...
    POL,FLT,WND,AT,clPick.residual_error(:),LUD,'VariableNames',varNames);

% NOTE: Can append new column via the following:
% ATT.newVariable = [];
