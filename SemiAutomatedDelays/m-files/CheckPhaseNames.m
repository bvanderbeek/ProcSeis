%% Check Phase Names
close all;
clear;
clc;

% Input
dataDir      = '~/research/CASCADIA/Waveforms_S'; % Data directory
theEvent     = [dataDir,'/event_cascadia.mat']; % Event structure
theStation   = [dataDir,'/station_cascadia.mat']; % Station structure
theATT       = [dataDir,'/ATT_cascadia.mat']; % The arrival time table
refModel     = 'ak135';

% Load data
load(theEvent,'event');
load(theStation,'station');
load(theATT,'ATT');

% Identify initial picks
iog = ATT.window == 0;
iog = true(size(iog)); % Update all!
N   = sum(iog);

% Get measured arrival times
TO = datenum(ATT.originTime(iog));
TA = datenum(ATT.arrivalTime(iog));
TT = (TA - TO)*24*60*60;

% Create source-reciever distance and event depth arrays
[~,ista] = ismember(ATT.station(iog),station.name);
[~,ievt] = ismember(ATT.event(iog),event.id);
DLT = distance(event.latitude(ievt),event.longitude(ievt),...
    station.latitude(ista),station.longitude(ista));
DPT = event.depth(ievt);
PHS = ATT.phase(iog);

% Loop over arrivals
nbad = 0; % Number of invalid phase names
nmul = 0; % Number of multiple arrivals
tic
for ii = 1:N
    % Predicted travel-time for defined phase
    [tt1D,~,~,~,outPhase,iwarn] = taup_time(refModel,PHS{ii},DLT(ii),DPT(ii),0);
    
    % Check for multiples
    if iwarn == 1
        nbad = nbad + 1;
        PHS{ii} = outPhase;
    elseif iwarn == 2
        nmul = nmul + 1;
    end
    
    % Valid arrival returned?
    if isinf(tt1D) || isnan(tt1D)
        error('Unexpected behavior');
        % Define phases to check
        if strcmp(PHS{ii},'S')
            phs_list = {'Sdiff'}; %,'SKS','SKKS'};
        elseif strcmp(PHS{ii},'sS')
            phs_list = {'sSdiff'}; %,'sSKS','sSKKS'};
        else
            error('Unknown phase.');
        end
        
        % Check alternative phases
        jmin = 0;
        dt   = 1e6;
        for jj = 1:length(phs_list)
            tt1D = taup_time(refModel,phs_list{jj},DLT(ii),DPT(ii),0);
            if abs(TT(ii) - tt1D) < dt
                jmin = jj;
                dt   = abs(TT(ii) - tt1D);
            end
        end
        
        % Update phase if a valid arrival was found
        if jmin == 0
            error('No valid arrival found');
        else
            PHS{ii} = phs_list{jmin};
        end
    end
end
toc

% Update phase names
% ATT.phase(iog) = PHS;