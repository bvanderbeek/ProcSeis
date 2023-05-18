%% Make Priority Event List
% Arranges events by year and number of arrivals
close all;
clear;
clc;

% Input
dataDir      = '~/research/CASCADIA/Waveforms_S'; % Data directory
theEvent     = [dataDir,'/event_cascadia.mat']; % Event structure
theATT       = [dataDir,'/ATT_cascadia.mat']; % The arrival time table

% Load data structures
load(theEvent,'event');
load(theATT,'ATT');

% Identify original arrivals
iog = ATT.window == 0;

% Compute number of arrivals for each event
EID      = unique(ATT.event(iog));
[~,ievt] = ismember(ATT.event(iog),EID);
NEVT     = accumarray(ievt,1,size(EID),@sum,0);
% Get event years
[~,ievt] = ismember(EID,event.id);
EYR      = event.originDate.Year(ievt);
% Unique list of years
YR = unique(event.originDate.Year);
YR = YR(:)';

% Initialize Event Queue
P  = [EID,EYR,NEVT];
Q  = zeros(size(P));
iq = 0;
while ~isempty(P)
    for yyyy = YR
        % Select events for specific year
        iin = P(:,2) == yyyy;
        nmax = max(P(iin,3));
        if ~isempty(nmax)
            iadd = find(iin & (P(:,3) == nmax),1,'first');
            if ~(length(iadd) == 1)
                error('Unexpected behavior');
            end
            iq   = iq + 1;
            Q(iq,:)   = P(iadd,:);
            P(iadd,:) = [];
        end
    end
end

%% Write priority queue
fid = fopen('EventPriorityQueue.txt','w');
for n = 1:size(Q,1)
    fprintf(fid,'%6d %6d %6d \n',Q(n,:));
end
fclose(fid);
