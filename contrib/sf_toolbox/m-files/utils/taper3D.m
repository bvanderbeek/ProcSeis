function K = taper3D(nx,ny,nz,numx,numy,numz,tf_hard,fbuff)
% Creates tappering array.
% Taper in third dimension is only at top
% Apply taper by G = F + K*(Fedge - F)

% Taper X
nb  = round(fbuff*numx);
tap = (1 + cos(2*pi*linspace(0,2*(numx-nb),2*(numx-nb)+1)./(2*(numx-nb))))./2;
tap = cat(2,ones(1,nb),tap,ones(1,nb));
% Hard taper?
if tf_hard
    tap(tap > 0) = 1;
end
Kx = [tap(1:numx)';zeros(nx-2*numx,1);tap(numx+2:end)'];

% Taper Y
nb  = round(fbuff*numx);
tap = (1 + cos(2*pi*linspace(0,2*(numy-nb),2*(numy-nb)+1)./(2*(numy-nb))))./2;
tap = cat(2,ones(1,nb),tap,ones(1,nb));
% Hard taper?
if tf_hard
    tap(tap > 0) = 1;
end
Ky = [tap(1:numy),zeros(1,ny-2*numy),tap(numy+2:end)];

% Taper Z
nb  = round(fbuff*numx);
tap = (1 + cos(2*pi*linspace(0,2*(numz-nb),2*(numz-nb)+1)./(2*(numz-nb))))./2;
tap = cat(2,ones(1,nb),tap,ones(1,nb));
% Hard taper?
if tf_hard
    tap(tap > 0) = 1;
end
Kz = [tap(1:numz)';zeros(nz-numz,1)];

% Combine
Kxy = repmat(sqrt((1-Kx)*(1-Ky)),1,1,nz);
Kz  = permute(repmat(1-Kz,1,nx,ny),[2,3,1]);
K   = 1-sqrt(Kxy.*Kz.*Kxy);

