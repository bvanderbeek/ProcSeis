function w = get_window(nsamp,fs,ttaper,varargin)
% GET_WINDOW: Creates a simple time domain window.
%
% INPUT
%     nsamp: Number of samples in signal
%        fs: Sampling frequency of signal
%    ttaper: Length of hanning taper in units of 1/fs at edge of window.
%            Use 0 for a boxcar.
%
% <optional>
%    tstart: Start time of signal
%    wstart: Time of start of window
%   wlength: Length of window
%
% OUTPUT
%         w: time domain window (nsamp x 1 array)
%
% NOTES
% + Without variable input arguments, just tapers ends of signal
%

% Check variable input arguments
if isempty(varargin)
    tstart  = [];
    wstart  = [];
    wlength = [];
elseif length(varargin) == 3
    tstart  = varargin{1};
    wstart  = varargin{2};
    wlength = varargin{3};
else
    error('Incorrect number of variable input arguments')
end

% Check inputs
if ~isempty(tstart)
    if ceil((wlength + 2*ttaper)*fs) > nsamp
        error('Length of window and taper exceeds signal length.');
    end
    if (wstart-ttaper) < tstart
        error('Requested window begins before start of trace');
    end
else
    if ceil(ttaper*fs) > nsamp
        error('Requested taper exceeds signal length.')
    end
end

if isempty(tstart)
    % Just tapering ends of trace to zero
    ntap = 1 + round(ttaper*fs);
    w    = hann(2*ntap);
    w    = [w(1:ntap); ones(nsamp-2*ntap,1); w((ntap+1):end)];
else
    % Creating window
    nwin   = 1 + round(wlength*fs);
    ntap   = 1 + round(ttaper*fs);
    nstart = 1 + round((wstart - tstart)*fs) - ntap;
    w      = hann(2*ntap);
    w      = [zeros(nstart,1); w(1:ntap); ones(nwin,1); w((ntap+1):end)];
    nadd   = nsamp-nstart-nwin-2*ntap;
    if nadd < 0
        warning('End of window exceeds signal length.');
        w = w(1:nsamp);
    else
        w = cat(1,w,zeros(nadd,1));
    end
end
