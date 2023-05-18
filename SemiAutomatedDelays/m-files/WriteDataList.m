%% Write Data Lists
% This script does the following:
% 1) Creates directories for every event with picks
% 2) Inside each event directory, a list of stations recording the event
%    followed by estimated arrival time is stored. This list is used to
%    retrieve waveform data via a Python script.
% Should only have to run this script once...
close all;
clear;
clc;

% Where to create event directories
dataDir = '~/research/CASCADIA/TEST'; %Waveforms_S';
% Arrival time table
theATT = '../data/ATT_cascadia.mat';

%% Create Directories and Arrival Lists
load(theATT,'ATT');

% Get unique event list
EID  = unique(ATT.event);
nevt = length(EID);

% Write data
for ii = 1:nevt
    % Make event directory
    evtDir = [dataDir,'/EVT',num2str(EID(ii))];
    mkdir(evtDir);
    
    % Get list of stations recording this event
    nrow = find(ATT.event == EID(ii));
    
    % AD-HOC FIX: Only event ID 697 has multiple phases recorded (S and
    % sS). The most observations and smallest errors are for S and so we
    % will only consider these at this time.
    if EID(ii) == 697
        nrow = nrow(strcmp(ATT.phase(nrow),'S'));
    end
    
    % Check unique station list returned
    if ~(length(unique(ATT.station(nrow))) == length(nrow))
        error('Non-unique stations!');
    end
    nrow = nrow(:)'; % Column vector for loop indices
    
    % Create arrival list
    fid = fopen([evtDir,'/ArrivalList.dat'],'w');
    for jj = nrow(1:end-1)
        fprintf(fid,'%s %s %s %s %s',ATT.network{jj},ATT.station{jj},...
            ATT.arrivalTime(jj),ATT.phase{jj});
        fprintf(fid,'\n');
    end
    % Last entry
    jj = nrow(end);
    fprintf(fid,'%s %s %s %s %s',ATT.network{jj},ATT.station{jj},ATT.arrivalTime(jj),ATT.phase{jj});
    fclose(fid);
end
