function [delays,tpick,ibad,CCF,itol,S,seis] = ICCS(tstart,fs,seis,ttaper,wstart,wlength,ccf_min,PickParams,maxlag,tshift,varargin)
% ICCS: Aligns seismic traces by iteratively cross-correlating individual
% traces with an array stack to compute time shifts. A simplified
% implementation of method described by Lou et al. (2013) for teleseismic
% delay time measurement.
%
% INPUT
%       tstart: Start time of seismic traces
%           fs: Sampling frequency of traces
%         seis: Array of seismic traces to be processed
%              - Dimensions are n-receivers by m-samples
%               - Assumed that traces have already been aligned on a 1D
%                 prediction
%       ttaper: Length of taper aadded to windowed seismograms (can be 0 for boxcar)
%       wstart: Start time of cross-correlation window
%      wlength: Length of cross-correlation window
%      ccf_min: Reject traces with correlation coefficient less than this
%               value from the stack
%   PickParams: Structure defining parameters for picking phase onset.
%               There are three options:
%               PickParams.method = 'manual'
%               + Manually pick array stack; no additional structure fields
%                 required
%               PickParams.method = 'stalta'
%               + Identify onset time from max of short-term over long-term
%                 average. Requires these additional fields,
%                 .
%               PickParams.method = 'lomax'
%               + Identify onset time using the filter picker designed by
%               Lomax et al. (2012). Requires these additional fields,
%                 .
%               PickParams.method = 'nopick'
%               + Use a constant user-defined onset time. Requires these 
%                 additional fields,
%                 .tpick: Onset time of arrival
%
% <optional>
%    filtopts: Vector describing Butter worth filter options defined as:
%              [min. frequency, max frequency, order, tf_zerophase]
%              - Default is no filter
%   niter_max: Maximum number of iterations to align traces
%              - Default is 10
%     tol_max: Tolerance for convergence
%              - Default is 2x the sampling frequency
%
% OUTPUT
%   delays: Trace delay times with respect to stacked trace
%      CCF: The maximum cross-correlation coefficient for each trace at the
%           final iteration
%     itol: Maximum delay calculated at each iteration
%        S: Stacked trace at final iteration
%     seis: Filtered (if applicable) and time-shifted seismograms
%
% NOTES
% Trace amplitudes are normalized by RMS amplitude in window before
% cross-correlating.
%
% REFERENCES
% [1] Lou, X., van der Lee, S., & Lloyd, S. (2013). AIMBAT: A python/matplotlib 
%     tool for measuring teleseismic arrival times. Seismological Research 
%     Letters, 84(1), 85-93.


% Assign options input arguments
if isempty(varargin)
    filtopts  = [];
    niter_max = 10;
    tol_max   = 2/fs;
elseif length(varargin) == 1
    filtopts = varargin{1};
elseif length(varargin) == 2
    filtopts = varargin{1};
    tol_max  = varargin{2};
elseif length(varargin) == 3
    filtopts  = varargin{1};
    tol_max   = varargin{2};
    niter_max = varargin{3};
else
    error('Incorrect number of optional input arguments.')
end

% Hidden option to run a test
tf_test = false;
if strcmp(tstart,'test')
    tf_test = true;
    % Synthetic parameters
    dt_max = 2; % Maximum magnitude of time shift added to trace
    N      = 10; % Number of synthetic traces to create
    Tc     = 10; % Central frequency of a Ricker wavelet
    Ts     = Tc/100; % Sampling period
    t      = -2*Tc:Ts:Tc*2; % Time vector
    
    % Build seismograms (Ricker Wavelets)
    seis   = zeros(N,length(t));
    dt_add = 2*dt_max*(rand(N) - 0.5);
    for ii = 1:N
        seis(ii,:) = (1 - 2*(pi^2)*(1/(Tc^2))*((t+dt_add(ii)).^2)).*exp(-(pi^2)*(1/(Tc^2))*((t+dt_add(ii)).^2));
    end
    
    % Define inputs
    tstart   = t(1);
    fs       = 1/Ts;
    ttaper   = 2;
    wstart   = -(Tc + dt_max);
    wlength  = 2*(Tc + dt_max);
    filtopts = [];
    tol_max  = 2/fs;
    
    % Make initial plot of seismograms
    figure(101);
    clf; hold on;
    for ii = 1:N
        plot(t,ii + seis(ii,:),'-k');
    end
    box on; grid on;
end

% Make sure we can find the toolbox for Lomax filter picker
if strcmp(PickParams.method,'lomax')
    if isempty(getenv('LOMAX'))
        warning('Missing LOMAX environment variable that points to filter picker toolbox. Defaulting to manual pick identification.');
        clear PickParams;
        PickParams.method = 'manual';
    else
        addpath(getenv('LOMAX'));
    end
end

% Identify number of traces and samples
[ntrace,nsamp] = size(seis);

% Define maximum lag
if isempty(maxlag)
    maxlag = nsamp - 1;
else
    maxlag = maxlag*fs;
end

% Time vector
t = tstart + linspace(0,nsamp-1,nsamp)./fs;

% Pre-processing
% Remove an initial delay if provided and get trace amplitudes
a = ones(ntrace,1);
if ~isempty(tshift) && (sum(abs(tshift)) > 10*eps)
    for ii = 1:ntrace
        seis(ii,:) = interp1(t - tshift(ii),seis(ii,:),t,'linear',0);
        w = get_window(nsamp,fs,ttaper,tstart,tshift(ii) + wstart,wlength);
        a(ii) = rms(seis(ii,w > 0));
    end
end
% Taper ends of traces
w    = get_window(nsamp,fs,ttaper);
seis = seis.*repmat(w(:)',ntrace,1);
% Filter traces
if ~isempty(filtopts)
    for ii = 1:ntrace
        seis(ii,:) = ButterFilter(seis(ii,:),filtopts(1:2),filtopts(3),fs,filtopts(4));
    end
end

% Initialize while-loop parameters
tol    = 100*tol_max;
niter  = 0;
itol   = [];
delays = zeros(ntrace,1);

% Run ICCS
ibad = false(ntrace,1); % Initially include all traces
% Add loop over window lengths?
while (tol > tol_max) && (niter < niter_max)
    % Update counter
    niter = niter + 1;
    
    % Stacked trace
    % Should I normalize trace amplitudes before stacking? For my synthetic
    % teleseisms, amplitudes tend not to vary significantly across array so
    % probably okay to skip for now.
    S = sum(seis(~ibad,:)./repmat(a(~ibad),1,nsamp),1)./sum(~ibad);
    % S = sum(seis(~ibad,:),1)./sum(~ibad); % Stacked trace
    
    % Pick stacked trace
    if strcmp(PickParams.method,'lomax')
        % Lomax Filter Picker Method
        % Reduce tolerance until pick is made
        tpick = [];
        dE    = 0;
        E1    = PickParams.E1;
        E2    = PickParams.E2;
        while isempty(tpick) && (E1 > 0) && (E2 > 0)
            % Filter picker parameters
            stw   = round(fs*PickParams.stw);
            ltw   = round(fs*PickParams.ltw);
            ptw   = round(fs*PickParams.ptw);
            E1    = PickParams.E1-dE;
            E2    = PickParams.E2-dE;
            % Make pick
            tpick = filterPicker(1/fs,S(:)',[stw ltw ptw E1 E2],[1 nsamp]);
            tpick = t(1) + tpick;
            % Update pick tolerances
            dE    = dE + 1;
        end
        if (dE - 1) > 0
            warning(['Tolerance reduced by ',num2str(dE-1),' before a pick was identified.'])
        end
        if isempty(tpick)
            warning('No pick identified. Setting arrival to 0 s.');
            tpick  = 0;
        elseif length(tpick) > 1
            % This is an ad-hoc selection in the case that multiple picks are
            % found
            warning('Multiple phases identified. Chosing earliest arrival...')
            tpick  = tpick(1);
        end
    elseif strcmp(PickParams.method,'stalta')
        % STA/LTA Method
        tpick = sta_lta_pick(S,fs,t(1),PickParams.stw,PickParams.ltw,PickParams.ttaper,PickParams.pick_window);
    elseif strcmp(PickParams.method,'manual')
        figure(101); clf;
        plot(t,S,'-b','linewidth',2);
        grid on;
        disp('Please make a pick...');
        shg;
        [tpick,~] = ginput(1);
        close(figure(101));
    elseif strcmp(PickParams.method,'nopick')
        % Use a constant pick time
        tpick = PickParams.tpick;
    end
    
    % Update arrival window
    w = get_window(nsamp,fs,ttaper,tstart,tpick + wstart,wlength);
    w = w(:)';
    
    % Normalize final stacked trace amplitude
    S  = S./rms(S(w > 0));
    
    % Calculate time shifts
    CCF = zeros(ntrace,1);
    DT  = zeros(ntrace,1);
    for ii = 1:ntrace
        % Cross-correlate individual traces with stacked-windowed trace
        si       = seis(ii,:);
        a(ii)    = rms(si(w > 0)); % Update rms amplitude in window
        si       = si./a(ii); % Scale trace
        [r,lags] = xcorr(si.*w,S.*w,maxlag); % Window both traces and xcorr
        lags     = lags./fs;
        keep     = (r == max(r));
        DT(ii)   = mean(lags(keep));  % A crude choice for delay?
        CCF(ii)  = max(r);
        
        % Remove delay from the traces via interpolation
        seis(ii,:) = interp1(t - DT(ii),seis(ii,:),t,'linear',0);
    end
    % Scale CCF by auto-correlation
    % CCF = CCF./max(xcorr(S.*w,S.*w));
    CCF = CCF./sum((w.*S).^2);
    
    % Identify bad traces
    ibad = CCF < ccf_min;
    if sum(ibad) == ntrace
        warning('No traces meet quality criteria.');
        niter = niter_max;
    end
    
    if any(abs(DT) > (maxlag/fs))
        warning('Measured delays reached maximum lag allowed.');
    end
    
    % Total time lags
    delays = delays + DT;
    
    % Update tolerance
    tol  = max(abs(DT(~ibad))); % Only include good traces for convergence criteria
    itol = cat(1,itol,tol);
end

% Display result
if niter < niter_max
    disp(['Converged at iteration ',num2str(niter),'.']);
else
    disp(['Reached iteration ',num2str(niter),' without converging.']);
end

% Plot test results
if tf_test
    figure(201);
    clf; hold on;
    for ii = 1:N
        plot(t,ii + seis(ii,:),'-k');
    end
    plot(t,ii + 1 + S,'-b','linewidth',2);
    box on; grid on;
end
