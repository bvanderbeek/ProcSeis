function PAZ = calc_polarization(S1,S2,S3)
% Returns polarization azimuth measured clockwise from channel 2
N = length(S1);
% Covariance matrix of stack
c11   = sum((S1 - mean(S1)).^2)./N;
c12   = sum((S1 - mean(S1)).*(S2 - mean(S2)))./N;
c13   = sum((S1 - mean(S1)).*(S3 - mean(S3)))./N;
c22   = sum((S2 - mean(S2)).^2)./N;
c23   = sum((S2 - mean(S2)).*(S3 - mean(S3)))./N;
c33   = sum((S3 - mean(S3)).^2)./N;

% Polarization analysis
[V,D] = eig([c11 c12 c13; c12 c22 c23; c13 c23 c33]);
D     = diag(D);
imax  = D == max(D);
V     = V(:,imax);
% V     = sign(V(2))*V; % Force upper quadrant
PAZ   = atan2d(V(1),V(2));
