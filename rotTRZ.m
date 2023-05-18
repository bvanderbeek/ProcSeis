function [seis,Rmat] = rotTRZ(seis,BAZ)

% Dimensions of seismic data array
[ntrace,nchan,~] = size(seis);

if (nchan > 1) && (nchan < 4)
    % Apply rotations to each trace
    Rmat = zeros(ntrace,nchan*nchan);
    for ii = 1:ntrace
        % Define rotation matrix. Note that back-azimuth is measured
        % clockwise from north but rotation matrix is defined such that
        % positive angles are counterclockwise.
        R = [cosd(BAZ(ii)),-sind(BAZ(ii)), 0;...
            sind(BAZ(ii)),  cosd(BAZ(ii)), 0;...
            0,              0            , 1];
        % Subset for number of channels
        R = R(1:nchan,1:nchan);
        % Apply rotation
        seis(ii,:,:) = R*squeeze(seis(ii,:,:));
        % Store rotation matrix components
        Rmat(ii,:) = R(:)';
    end
else
    error(['Cannot rotate to ''TRZ'' with ',num2str(nchan),' channels.']);
end
