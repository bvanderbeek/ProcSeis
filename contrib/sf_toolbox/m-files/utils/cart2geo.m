function [LON,LAT,R] = cart2geo(X,Y,Z,lon0,lat0,tf_elev)

% Typically, Z will be given as depth below Earth's surface (i.e.
% elevation).

% Earth Radius
Re = 6371;

% Get local z-coordinate
if tf_elev
    % If Z is elevation with respect to Earth's surface
    Z = sqrt(((Re+Z).^2) - (X.^2) - (Y.^2));
else
    % If Z is depth (negative) in cartesian box with surface at Z = 0
    % NOTE!!! This was the only Z definition in the old version
    Z = Re + Z;
end

% Store input array sizes
[ni,nj,nk] = size(Z);

% Convert cartesian box coordinates to latitude/longitude
% If we place the Earth in a box, the X,Y-plane falls along the equator. In
% this global cartesian system, the local SPECFEM x, y, and z coordinates
% for a model centered on the equator/prime meridian correspond to ygloal,
% zglobal, and xglobal respectively. In this geometry latitude and
% longitude can be obtained from a cartesian to spherical coordinate
% transform. Assumes perfectly spherical Earth.

% Get global cartesian coordinates
xyz_global = rotz(lon0)*roty(-lat0)*[Z(:)';X(:)';Y(:)'];

% Convert to geographic coordinates
[LON,LAT,R] = cart2sph(xyz_global(1,:)',xyz_global(2,:)',xyz_global(3,:)');
LON = rad2deg(LON);
LAT = rad2deg(LAT);

% Return arrays in same dimension
LON = reshape(LON,ni,nj,nk);
LAT = reshape(LAT,ni,nj,nk);
R   = reshape(R,ni,nj,nk);
