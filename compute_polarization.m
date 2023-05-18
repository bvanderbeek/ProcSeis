function [azm,elv,D,V] = compute_polarization(u1,u2,u3)

% Demean each component
u1 = u1 - mean(u1);
u2 = u2 - mean(u2);
u3 = u3 - mean(u3);

% Covariance matrix components
N   = length(u1);
c11 = sum(u1.*u1)/N;
c12 = sum(u1.*u2)/N;
c13 = sum(u1.*u3)/N;
c22 = sum(u2.*u2)/N;
c23 = sum(u2.*u3)/N;
c33 = sum(u3.*u3)/N;

% Solve eigen-value problem for polarization
[V,D] = eig([c11 c12 c13; c12 c22 c23; c13 c23 c33],'vector');

% Sort such that channel 3 has most energy
[D,iD] = sort(abs(D));
V      = V(:,iD);

% Compute polarization angles
[azm,elv] = cart2sph(V(1,3),V(2,3),V(3,3));
