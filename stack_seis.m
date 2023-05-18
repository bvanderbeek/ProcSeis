function [S,seis] = stack_seis(t,seis,DT)

% Dimensions of seismic data array
[ntrace,nsamp] = size(seis);

% Coordinate vectors
t      = t(:)';
itrace = (1:ntrace)';

% Interpolate seismograms to common phase-aligned time vector
seis = interp2(t,itrace,seis,repmat(t,ntrace,1) + repmat(DT(:),1,nsamp),...
    repmat(itrace,1,nsamp),'linear',0);

% Stack
S = mean(seis,1);