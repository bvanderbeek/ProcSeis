function [delays,S,ddt,dtij,A,res] = MCCS(tstart,fs,seis,dtpick,ttaper,wstart,wlength,maxlag,varargin)

% Identify number of traces and samples
[ntrace,nsamp] = size(seis);

% Define maximum lag
if isempty(maxlag)
    maxlag = nsamp - 1;
else
    maxlag = maxlag*fs;
end

% Define option input arguments
if isempty(varargin)
    wts       = ones(ntrace,1);
    filtopts  = [];
elseif length(varargin) == 1
    wts      = varargin{1};
    filtopts = [];
elseif length(varargin) == 2
    wts      = varargin{1};
    filtopts = varargin{2};
else
    error('Incorrect number of optional input arguments.')
end

% Pre-processing
% Taper ends of traces
w    = get_window(nsamp,fs,ttaper);
seis = seis.*repmat(w(:)',ntrace,1);
% Filter traces
if ~isempty(filtopts)
    for ii = 1:ntrace
        seis(ii,:) = ButterFilter(seis(ii,:),filtopts(1:2),filtopts(3),fs,filtopts(4));
    end
end


% Multi-channel cross-correlation
N    = sum(1:(ntrace-1));
row  = zeros(N,2); % Coefficient matrix row-index
col  = zeros(N,2); % Coefficient matrix column-index
aij  = zeros(N,2); % TO ADD: Coefficient matrix value weight
dtij = zeros(N,1); % Delay between the ith and jth stations
kk   = 0; %
for ii = 1:(ntrace-1)
    % Arrival window for ith-trace
    wi = get_window(nsamp,fs,ttaper,tstart,dtpick(ii)+wstart,wlength);
    wi = wi(:)';
    
    % Normalized ith-trace
    si = seis(ii,:);
    si = si./rms(si(wi > 0));
    
    for jj = (ii+1):ntrace
        % Update counter
        kk = kk + 1;
        
        % Arrival window for jth-trace
        wj = get_window(nsamp,fs,ttaper,tstart,dtpick(jj)+wstart,wlength);
        wj = wj(:)';
        
        % Normalized jth-trace
        sj = seis(jj,:);
        sj = sj./rms(sj(wj > 0));
        
        % Cross-correlate windowed traces
        [r,lags] = xcorr(si.*wi,sj.*wj,maxlag);
        keep     = (r == max(r)); % THIS IS PERHAPS A CRUDE ESTIMATE
        
        % How to define weight?
        %wij = wts(ii); % Weight by the trace currently being cross-correlated with the remaining?
        %wij = sqrt(wts(ii)*wts(jj)); % Compsit weight option 1?
        wij = sqrt((wts(ii)^2) + (wts(jj)^2)); % Composit weight option 2?
        
        % Weighted delay time Ti - Tj
        dtij(kk) = wij*mean(lags(keep))/fs; % THIS IS PERHAPS A CRUDE ESTIMATE
        
        % Weighted coefficient
        aij(kk,:) = wij*[1 -1];
        
        % Indexing
        row(kk,1) = kk;
        row(kk,2) = kk;
        col(kk,1) = ii;
        col(kk,2) = jj;
        
    end
    %disp(ii);
end
% Coefficient matrix
A = sparse(row(:),col(:),aij(:),N,ntrace);
% Add weighted zero mean constraint
A    = cat(1,A,wts(:)');
dtij = cat(1,dtij,0);

% Best-fit delay times
delays = lsqr(A,dtij,[],1000);
res    = dtij(:) - A*delays(:); % Residual vector for error estimation

% Create a final stack
t = tstart + linspace(0,nsamp-1,nsamp)./fs;
for ii = 1:ntrace
    % Remove delay from the trace via interpolation
    seis(ii,:) = interp1(t - delays(ii),seis(ii,:),t,'linear',0);
end
w   = (t >= 0) & (t <= wlength);
amp = rms(seis(:,w),2);
S   = sum(seis./repmat(amp,1,nsamp),1)./ntrace; % Amplitude-scaled stack
% S = sum(seis,1)./ntrace; % Unscaled stack

% Estimate errors using Eq. 8 of Vandecar and Crosson, 1990. This is an underestimate!
ddt  = zeros(ntrace,1);
for ii = 1:length(ddt)
    % 1) Error Estimate following Vandecar and Crosson (1990)
    rji = 0;
    rij = 0;
    % Sum of squared residuals associated with ith-trace
    for ji = 1:(ii-1)
        ind = sum((ntrace-(1:ii))) + (ji-ntrace);
        rji = rji + (res(ind)^2);
    end
    for ij = (ii+1):ntrace
        ind = sum((ntrace-(1:ii))) + (ij-ntrace);
        rij = rij + (res(ind)^2);
    end
    % Estimated standard error
    ddt(ii) = sqrt((1/(ntrace-2))*(rji + rij));
end

if any(abs(delays) > (maxlag/fs))
    warning('Measured delays reached maximum lag allowed.');
end
