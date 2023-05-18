function [A,lwin] = get_rms_amp(t,seis,wstart,wlength)

% Get number of traces
[ntrace,nsamp] = size(seis);

% Vector of window lengths
if length(wlength) == 1
    wlength = wlength(ones(ntrace,1));
end
if ~(length(wlength) == ntrace)
    error('Input ''wlength'' must be a scalar or a vector the size of ''wstart''.');
end

% Get RMS trace amplitudes in window
A    = ones(ntrace,1);
lwin = true(ntrace,nsamp);
for ii = 1:ntrace
    ll         = (t >= wstart(ii)) & (t <= (wstart(ii) + wlength(ii)));
    A(ii)      = rms(seis(ii,ll));
    lwin(ii,:) = ll;
end
