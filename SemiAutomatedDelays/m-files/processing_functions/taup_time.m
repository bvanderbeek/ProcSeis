function [tt,rayP,inc,toff,outPhase,iwarn] = taup_time(aModel,aPhase,delta,sdepth,rdepth)
% TAUP_PATH: A very simple matlab wrapper to call the taup_time function
% from the TauP Java toolbox.
%
% INPUT
%     aModel: A TauP recognized model name, e.g. 'iasp91' (string)
%     aPhase: A TauP recognized seismic phase name (string)
%      delta: A distance (deg.)
%     sdepth: Source depth (positive into Earth; km)
%     rdepth: Receiver depth (positive into Earth; km)
%
% OUTPUT
%         tt: Travel-time (s)
%       rayP: Ray parameter (s/deg)
%        inc: Ray incidence angle at station (deg. from vertical)
%       toff: Ray take-off angle at source (deg. from vertical)
%   outPhase: Output phase type. Should match aPhase. However, if aPhase is
%             'P' or 'S', then the first arriving compressional or shear
%             phase is returned.
%
% NOTES
% + The Java matlab object can compute a variety of ray path attributes
%   that may be added to this function if needed (e.g. pierce points).
%
% + I modified the behavior of this function when TauP predictes no
%   arrivals for requested phase. Will throw a warning and return infinte 
%   values. BPV JAN 2023.
%
% REFERENCES
% + Crotwell, H. P., T. J. Owens, and J. Ritsema (1999). The TauP ToolKit: 
%   Flexible Seismic Travel-Time and Raypath Utilities, Seismological 
%   Research Letters
%   http://www.seis.sc.edu
%
% B. VanderBeek AUG-2019
%


% Check if we need to add TauP jar. We don't want to always add it because
% it can be slow to add to the java class path.
jpath = javaclasspath('-dynamic');
if ~any(strcmp(jpath,getenv('TAUPJAR')))
    if isempty(getenv('TAUPJAR'))
        error('Missing ''TAUPJAR'' environment variable pointing to TauP jar file.');
    else
        javaaddpath(getenv('TAUPJAR'));
    end
end

% Initialize warning. (0) Consistent arrival found, (1) alternative arrival
% found, (2) multiple arrivals found, and (3) no arrival found.
iwarn = 0;

% Create the travel-time java object. Takes ~1 ms to make.
ttobj = javaObject('edu.sc.seis.TauP.TauP_Time',aModel);

% Set the TauP_Time parameters
ttobj.setPhaseNames({aPhase});
ttobj.setSourceDepth(sdepth);
ttobj.setReceiverDepth(rdepth);

% Calculate arrival times
ttobj.calculate(delta);

% Extract results
N = ttobj.getNumArrivals;

% If no arrival returned, check these alternatives
if N == 0
    % Phases that are part of the same travel-time branch
    switch aPhase
        case 'P'
            alt_phs = 'Pdiff';
        case 'pP'
            alt_phs = 'pPdiff';
        case 'S'
            alt_phs = 'Sdiff';
        case 'sS'
            alt_phs = 'sSdiff';
        otherwise
            alt_phs = '';
    end
    
    % Check alternative phase
    if ~isempty(alt_phs)
        warning(['No ',aPhase,' arrivals predicted at delta = ',num2str(delta),' deg. Looking for ',alt_phs,' arrival...']);
        ttobj.setPhaseNames({alt_phs});
        ttobj.calculate(delta);
        N     = ttobj.getNumArrivals;
        iwarn = 1;
    end
end

% Evaluate TauP object
if N > 0
    % Arrival(s) returned
    if N > 1
        iwarn = 2;
        warning(['Multiple ',aPhase,' arrivals predicted at delta = ',num2str(delta),' deg.; using earliest arrival.']);
    end
    ttout    = ttobj.getArrival(0); % Indexing starts at 0
    outPhase = char(ttout.getName);
    tt       = ttout.getTime;
    rayP     = ttout.getRayParam*pi/180;
    inc      = ttout.getIncidentAngle;
    toff     = ttout.getTakeoffAngle;
else
    % No arrival returned
    error(['No ',aPhase,' arrivals predicted at delta = ',num2str(delta),' deg.']);
%     iwarn    = 3;
%     tt       = NaN;
%     rayP     = NaN;
%     inc      = NaN;
%     toff     = NaN;
%     outPhase = aPhase;  
end
