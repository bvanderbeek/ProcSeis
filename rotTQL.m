function [seis,Rmat] = rotTQL(seis,BAZ,INC)

% Dimensions of seismic data array
[ntrace,nchan,~] = size(seis);

if nchan == 3
    % Apply rotations to each trace
    Rmat = zeros(ntrace,nchan*nchan);
    for ii = 1:ntrace
        % Define rotation matrix. Note that back-azimuth is measured
        % clockwise from north but rotation matrix is defined such that
        % positive angles are counterclockwise.
        % Rotation about x-axis
        R = [1,              0,               0;...
             0, cosd(-INC(ii)), -sind(-INC(ii));...
             0, sind(-INC(ii)),  cosd(-INC(ii))];
        % Rotation about z-axis
        R = R*[cosd(BAZ(ii)), -sind(BAZ(ii)), 0;...
               sind(BAZ(ii)),  cosd(BAZ(ii)), 0;...
               0,              0,             1];
        % Apply rotation
        seis(ii,:,:) = R*squeeze(seis(ii,:,:));
        % Store rotation matrix
        Rmat(ii,:) = R(:)';
    end
else
    error('Rotation to ''TQL'' requires three channels.');
end
