function [DT,ddt,CCF] = SCC(t,fs,seis,dt_trace,tf_bad,wlength,tbuff,ttaper,tf_norm)
% SCC: Aligns seismic traces by cross-correlating individual traces with
% an array stack to compute time shifts. A simplified implementation of
% the method described by Lou et al. (2013).
%
% INPUT
%          t: Reduced time vector for seismograms (i.e. zero time is
%             assumed to correspond to the phase arrival time).
%	      fs: Sampling frequency of traces.
%       seis: Array of seismic traces to be processed.
%              - Dimensions are n-receivers by m-samples
%              - Assumes that traces have already been aligned on a 1D
%                prediction
%	dt_trace: A priori delays times for each trace such that t +
%	          dt_trace(ii) = 0 at onset of seismic phase of interest.
%     tf_bad: Logical vector of traces to exclude from analysis. Traces are
%            excluded when tf_bad is 'true'.
%    wlength: Length of phase window (1/fs) that starts at t = 0.
%      tbuff: A 2-element vector defining a buffer to add before and after 
%             the start/end of phase window (1/fs; always positive).
%     ttaper: Length of hanning tapper added to ends of phase window (1/fs).
%
% <Description of Window Parameters>
%   /|---|------|----|\
%  / |   |      |    | \
% /  |   |      |    |  \
%   -tb1 0      wlen wlen+tb2
% Phase window starts at 0 and extends to wlength. The window may be
% extended by tb1 (i.e. tbuff(1)) prior to t = 0 and by tb2 (i.e. tbuff(2))
% following the end of the phase window. Additionally, a symmetric tapper
% may be added to ends of phase window.
%
% OUTPUT
%       DT: Trace delay times with respect to reduced time vector
%      CCF: The maximum cross-correlation coefficient for each trace
%
% NOTES
%
%
% REFERENCES
% [1] Lou, X., van der Lee, S., & Lloyd, S. (2013). AIMBAT: A python/matplotlib 
%     tool for measuring teleseismic arrival times. Seismological Research 
%     Letters, 84(1), 85-93.


% Identify number of traces and samples
[ntrace,nsamp] = size(seis);

% Define maximum lag
maxlag = round((wlength + sum(tbuff))*fs);

% Interpolate traces to align phases on current delays
seis = seis';
seis = interp2(1:ntrace,t(:),seis,repmat(1:ntrace,nsamp,1),...
    repmat(t,1,ntrace) + repmat(dt_trace(:)',nsamp,1),'linear',0);

% Define stacked waveform
if tf_norm
    % Amplitude-normalized stack
    lwin = (t >= 0) & (t <= wlength);
    % Amplitude normalized traces
    seis = seis./repmat(rms(seis(lwin,:),1),nsamp,1);
    % Amplitude normalized stacked trace
    S    = sum(seis(:,~tf_bad),2)./sum(~tf_bad);
    S    = S./rms(S(lwin));
else
    % Straight stack
    S = sum(seis(:,~tf_bad),2)./sum(~tf_bad);
end
S    = S';
seis = seis';

% Phase Window
w = get_window(nsamp,fs,ttaper,t(1),-tbuff(1),wlength + sum(tbuff));
w = w(:)';

% Apply window
S    = w.*S;
seis = repmat(w,ntrace,1).*seis;

% Calculate time shifts
CCF = zeros(ntrace,1);
DT  = zeros(ntrace,1);
for ii = 1:ntrace
    if ~tf_bad(ii)
        % Normalized cross-correlation with stack
        [r,lags] = xcorr(seis(ii,:),S,maxlag,'coeff');
        % The time lags
        lags     = lags./fs;
        % Identify maximum correlation
        keep     = (r == max(r));
        % Define delay (average or select first maximum?)
        DT(ii)   = mean(lags(keep));
        CCF(ii)  = max(r);
    end
end

% Estimate errors based on similarity between traces and stacked waveform
% following Chevrot, GJI 2002.

% Re-define stacked trace aligning phases on new delays
seis = seis';
seis = interp2(1:ntrace,t(:),seis,repmat(1:ntrace,nsamp,1),...
    repmat(t(:),1,ntrace) + repmat(DT(:)',nsamp,1),'linear',0);
% Stack with optional normalization
if tf_norm
    lwin = (t >= 0) & (t <= wlength);
    % Amplitude normalized traces
    seis = seis./repmat(rms(seis(lwin,:),1),nsamp,1);
    % Amplitude normalized stacked trace
    S    = sum(seis(:,~tf_bad),2)./sum(~tf_bad);
    S    = S./rms(S(lwin));
else
    S = sum(seis(:,~tf_bad),2)./sum(~tf_bad);
end
S    = S';
seis = seis';

% Window seismograms
w    = get_window(nsamp,fs,ttaper,t(1),-tbuff(1),wlength + sum(tbuff));
S    = (w(:)').*S;
seis = repmat(w(:)',ntrace,1).*seis;

% Normalized auto-correlation of stack
AC = xcorr(S,S,'coeff');

ddt = zeros(ntrace,1);
for ii = 1:ntrace
    % Normalized cross-correlation with stack
    [r,lags] = xcorr(S,seis(ii,:),'coeff');
    
    % Identify where the auto-correlation coefficient exceeds the
    % correlation coefficient for the ith-trace. The width of this region
    % defines the error (note auto-correlation is symmetric so we only need
    % to consider one intersection).
    r  = max(r);
    jn = find((AC >= r),1,'last');
    
    % Estimated error
    ddt(ii) = (((lags(jn+1)-lags(jn))/(AC(jn+1)-AC(jn)))*(r-AC(jn)) + lags(jn))/fs;
end
