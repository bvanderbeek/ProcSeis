function S = align_traces(t,S,DT)
% Align traces via interpolation

% Seismogram array size
[nsta,nsamp,nchan] = size(S);

% Loop over channels
for ii = 1:nchan
    si = interp2(1:nsta,t(:),squeeze(S(:,:,ii))',repmat(1:nsta,nsamp,1),...
        repmat(t(:),1,nsta) + repmat(DT(:)',nsamp,1),'linear',0);
    S(:,:,ii) = si';
end