function [az,dt,L,tdelay,phi,L1,L2] = splitSC(U,fs,dt_max,nang,crit,tf_plot)
% SPLITRC: Measure shear-wave splitting via the Eigen-value method of
% Silver and Chan (1991).
%
% INPUT
%         U: Two-component signal (2xN array) 
%        fs: Sampling frequency of U (Hz)
%    dt_max: Maximum delay to consider (s)
%      nang: Number of angles in grid search on the interval [-pi/2, pi/2]
%      crit: Criteria for evaluating eigen-values ('L2/L1' or 'L1/L2')
%   tf_plot: If true, will make a plot of the eigen-value ratio surface
%
% OUTPUT
%       az: Fast azimuth(s) measured positive counter-clockwise from u1 (rad).
%       dt: Split time(s) (s)
%        L: Array of eigen-values ratios. The row index corresponds to each 
%           angle in the grid search while the column index corresponds to 
%           the delay
%   tdelay: Vector of delay times considered (s)
%      phi: Vector of angles considered in grid search (rad.)
%       L1: Maximum eigenvalue array
%       L2: Minimum eigenvalue array
%
% NOTES
%
%
% REFERENCES
% + Silver, P. G., & Chan, W. W. (1991). Shear wave splitting and 
%   subcontinental mantle deformation. Journal of Geophysical Research: 
%   Solid Earth, 96(B10), 16429-16454.
%
% B. VanderBeek (Nov-2020)
%

% Define maximum lag in samples
nsamp = size(U,2);
ndt   = 1 + floor(fs*dt_max);

% Define reference time vector
t = (linspace(0,nsamp-1,nsamp) - ((nsamp-1)/2))./fs;

% Grid search over fast directions and split times
phi    = linspace(-(pi/2),(pi/2),nang);
tdelay = linspace(0,ndt-1,ndt)/fs;
L1     = zeros(nang,ndt);
L2     = zeros(nang,ndt);
% Loop over potential fast directions
for iang = 1:nang
    % We are looking for the angle that places the fast-polarized wave onto
    % the x1-direction hence we use the transpose of the rotation matrix
    % (i.e. negative angle).
    R = [cos(phi(iang)), -sin(phi(iang)); sin(phi(iang)), cos(phi(iang))];
    R = R';
    W = R*U;
    
    % Loop over potential split times
    for jdt = 1:ndt
        % Remove delay via linear interpolation
        % Indexing
        tf   = t - (tdelay(jdt)/2);
        ts   = t + (tdelay(jdt)/2);
        indf = 1 + floor(fs*(tf - t(1)));
        indf = min(max(indf,1),nsamp-1);
        inds = 1 + floor(fs*(ts - t(1)));
        inds = min(inds,nsamp-1);
        
        % Interpolated signal with linear extrapolation. This is faster
        % than using Matlab's interpolation functions.
        f = (1 - fs*(tf - t(indf))).*W(1,indf) + fs*(tf - t(indf)).*W(1,indf+1);
        s = (1 - fs*(ts - t(inds))).*W(2,inds) + fs*(ts - t(inds)).*W(2,inds+1);
        % Constant extrapolation value
        iout    = (indf == 1);
        f(iout) = W(1,1);
        iout    = (inds == (nsamp-1));
        s(iout) = W(2,end);
        
        % Rotate back to starting coordinate system
        V = (R')*[f;s];
        
        % Compute covariance matrix
        c11 = sum((V(1,:) - mean(V(1,:))).^2)/nsamp;
        c22 = sum((V(2,:) - mean(V(2,:))).^2)/nsamp;
        c12 = sum((V(1,:) - mean(V(1,:))).*(V(2,:) - mean(V(2,:))))/nsamp;
        
        % Compute eigen-values used to assess linearity of particle motion.
        L1(iang,jdt) = (c11 + c22 + sqrt(((c11+c22).^2) - 4*((c11.*c22) - (c12.^2))))./2;
        L2(iang,jdt) = (c11 + c22 - sqrt(((c11+c22).^2) - 4*((c11.*c22) - (c12.^2))))./2;
    end
end

% Find orientation and lag based on linearity of unsplit waveform
switch crit
    case 'L2/L1'
        % Minimize L2/L1
        L          = (L2./L1);
        [iang,jdt] = find(L == min(L(:)));
    case 'L1/L2'
        % Maximize L1/L2
        L          = (L1./L2);
        [iang,jdt] = find(L == max(L(:)));
    case 'highbred'
        % Highbred
        L          = ((L2./L1) + (1 - (L1./max(L1(:)))))./2;
        [iang,jdt] = find(L == min(L(:)));
end

% Best fit parameters
az = phi(iang);
dt = tdelay(jdt);

% Plot linearity surface
if tf_plot
    figure; hold on;
    contourf(tdelay,phi*180./pi,log10(L),21);
    plot(dt,az*180./pi,'.r','markersize',25);
    box on; title(['Linearity Coefficient (log_1_0[',crit,'])']);
    xlabel('delay (s)');
    ylabel('azimuth (deg.)');
    axis square; colorbar;
end
