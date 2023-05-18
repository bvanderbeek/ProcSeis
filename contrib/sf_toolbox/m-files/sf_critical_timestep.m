function dt_crit = sf_critical_timestep(v,h,ngll,dim)
% Maximum time step scales positively with element size and inverse of wave speed
% This only provides an estimate. Always check mesh mesher output log for
% suggested time step.

% I copied this from utils/critical_timestep.m

% Tabulated values of critical frequency (non-dimensional, 1D)
% Omega_max(p) = Omega_max(ngll-1)
Omega_max = [2.0000000e+00 4.8989795e+00 8.6203822e+00 1.3540623e+01 1.9797952e+01 2.7378050e+01 ...
     3.6256848e+01 4.6421894e+01 5.7867306e+01 7.0590158e+01 8.4588883e+01 9.9862585e+01 ...
     1.1641072e+02 1.3423295e+02 1.5332903e+02 1.7369883e+02 1.9534221e+02 2.1825912e+02 2.4244948e+02];
 
% Stability factor for leapfrog timescheme
C = 2;

% Critical time step, assumes a cube element
dt_crit = C*h./v./sqrt(dim)/Omega_max(ngll-1);