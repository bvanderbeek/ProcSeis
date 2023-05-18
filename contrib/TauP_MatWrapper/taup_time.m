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
%   arrivals for requested phase. If, and only if, the requested phase is P
%   or S and no arrival is returned, check Pdiff and Sdiff. Otherwise,
%   error. BPV DEC 2022.
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

% Check if arrival was returned
if (N == 0) && (strcmp(aPhase,'P') || strcmp(aPhase,'S'))
    % Find first arrival that matches input phase
    warning(['No ',aPhase,' arrivals predicted at delta = ',num2str(delta),' deg. Looking for ',aPhase,'diff arrival...']);
    if strcmp(aPhase,'P')
        % Check Pdiff
        ttobj.setPhaseNames({'Pdiff'});
        ttobj.calculate(delta);
        N = ttobj.getNumArrivals;
    elseif strcmp(aPhase,'S')
        % Check Sdiff
        ttobj.setPhaseNames({'Sdiff'})
        ttobj.calculate(delta);
        N = ttobj.getNumArrivals;
    end
end

iwarn = 0;
if N > 0
    if N > 1
        iwarn = 2;
        warning(['Multiple ',aPhase,' arrivals predicted at delta = ',num2str(delta),' deg.; using first.']);
    end
    ttout    = ttobj.getArrival(0); % Indexing starts at 0
    outPhase = char(ttout.getName);
    tt       = ttout.getTime;
    rayP     = ttout.getRayParam*pi/180;
    inc      = ttout.getIncidentAngle;
    toff     = ttout.getTakeoffAngle;
else
    iwarn = 1;
    warning(['No ',aPhase,' arrivals predicted at delta = ',num2str(delta),' deg.']);
    tt       = Inf;
    rayP     = Inf;
    inc      = Inf;
    toff     = Inf;
    outPhase = aPhase;
end
