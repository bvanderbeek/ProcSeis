function h = ButterFilter(g,corners,order,fs,tf_zerophase)
% BUTTERFILTER--Simplifies action of applying a Butterworth (band-pass)
% filter to a signal. Also detrends the signal.
%
% INPUT
%              g: unfilterd signal
%        corners: corner frequencies of Butterworth filter [min,max] (Hz)
%          order: order of Butterworth filter
%             fs: sampling frequency (Hz)
%   tf_zerophase: If true, uses a zero-phase filter
%
% OUTPUT
%              h: filterd signal  
%

% Check dimensions
dim = size(g);
if (length(dim) ~= 2) && ~any(dim==1)
    error('Cannot do multidimensional arrays!');
end
% Force column vector
g = g(:);

% Check for valid corner frequencies
if any(2*corners > fs)
    warning('Corner frequency exceeds Nyquist Frequency! Data not filtered.')
    h = g;
else
    % Detrend and demean
    xvec = linspace(1,max(dim),max(dim))';
    P    = polyfit(xvec,g,1);
    g    = g - P(1)*xvec - P(2);
    g    = g - mean(g);
    
    % Taper
    ntap = ceil(fs*max(1./corners));
    tap  = hann(2*ntap);
    tap  = [tap(1:ntap); ones(length(g) - 2*ntap,1); tap((ntap+1):end)];
    g    = g.*tap;
    
    % Apply filter
    wlim  = sort(corners*2*(1/fs)); % use sort to idiot proof input
    [b,a] = butter(order,wlim);
    if tf_zerophase
        h = filtfilt(b,a,g);
    else
        h = filter(b,a,g);
    end
end

% Return h with same dimensions as g
if dim(1) == 1
    h = h';
end
