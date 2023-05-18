function [az,dt,CCF,tdelay,phi] = splitRC(U,fs,dt_max,nang,tf_plot)
% SPLITRC: Measure shear-wave splitting via rotation correlation method
% (i.e. Bowman & Ando, 1987).
%
% INPUT
%         U: Two-component signal (2xN array) 
%        fs: Sampling frequency of U (Hz)
%    dt_max: Maximum delay to consider (s)
%      nang: Number of angles in grid search on the interval [-pi/2, pi/2]
%   tf_plot: If true, will make a plot of the correlation surface
%
% OUTPUT
%       az: Fast azimuth(s) measured positive counter-clockwise from u1 (rad).
%       dt: Split time(s) (s)
%      CCF: Array of normalized cross-correlation values. The row index
%           corresponds to each angle in the grid search while the column
%           index corresponds to the delay
%   tdelay: Vector of delay times considered (s)
%      phi: Vector of angles considered in grid search (rad.)
%
% NOTES
%
%
% REFERENCES
% + Bowman, J. Roger, and Masataka Ando. "Shear-wave splitting in the 
%   upper-mantle wedge above the Tonga subduction zone." Geophysical 
%   Journal International 88.1 (1987): 25-41.
%
% B. VanderBeek (Nov-2020)
%

% Define maximum lag in samples
ndt = 1 + floor(fs*dt_max);

% Grid search over fast directions
phi    = linspace(-(pi/2),(pi/2),nang);
tdelay = linspace(-ndt,ndt,1+2*ndt)./fs;
CCF    = zeros(nang,1 + 2*ndt);
for iang = 1:nang
    % We are looking for the angle that places the fast-polarized wave onto
    % the x1-direction hence we use the transpose of the rotation matrix
    % (i.e. negative angle).
    R = [cos(phi(iang)), -sin(phi(iang)); sin(phi(iang)), cos(phi(iang))];
    R = R'; %%%%% GUARDA!
    W = R*U;
    
    % Compute normalized cross-correlation
    CCF(iang,:) = xcorr(W(1,:),W(2,:),ndt,'coeff');
end

% Identify orientation and lag based on maximum of cross-correlation.
[iang,jdt] = find(CCF == max(CCF(:)));
az         = phi(iang);
dt         = (jdt - (ndt + 1))./fs;

% If dt > 0, then fast and slow waveforms are negatively-correlated in FSZ
% system but positively-correlated in SFZ system. We always return fast
% direction but sense of correlation is indicated by sign of dt.
ianti      = dt > 0;
s          = sign(az);
s(s == 0)  = 1;
az(ianti)  = az(ianti) - s(ianti).*(pi/2);

% Plot correlation surface and waveforms
if tf_plot
    figure; hold on;
    contourf(tdelay,phi*180./pi,CCF,21);
    plot(-abs(dt),az*180./pi,'.r','markersize',25);
    box on; title('Correlation Coefficient');
    xlabel('delay (s)');
    ylabel('azimuth (deg.)');
    axis square; colorbar; caxis([-1 1]);
    
%     nsamp = size(U,2);
%     t     = (linspace(0,nsamp-1,nsamp) - (nsamp/2))./fs;
%     az0   = mean(az);
%     if mean(dt) > 0
%         az0 = (pi/2) + az0;
%     end
%     R     = [cos(az0), -sin(az0); sin(az0), cos(az0)];
%     R     = R'; %%%%% GUARDA!
%     Ufs   = R*U;
%     figure;
%     subplot(2,1,1); hold on;
%     plot(t,U(1,:),'-r','linewidth',2);
%     plot(t,U(2,:),'-b','linewidth',2);
%     box on; grid on; axis tight;
%     xlabel('time (s)');
%     legend('Radial','Transverse');
%     title('Transverse-Radial');
%     subplot(2,1,2); hold on;
%     plot(t,Ufs(1,:),'-b','linewidth',2);
%     plot(t,Ufs(2,:),'-r','linewidth',2);
%     box on; grid on; axis tight;
%     xlabel('time (s)');
%     legend('Fast','Slow');
%     title('Fast-Slow');
end
