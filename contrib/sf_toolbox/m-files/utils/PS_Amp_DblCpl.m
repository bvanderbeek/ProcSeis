function [up,ush,usv] = PS_Amp_DblCpl(strike,dip,rake,ray_phi,ray_theta)
% PS_AMP_DBLCPL: Compute the relative amplitude of P, Sh, and Sv waves
% generated from a double couple source mechanism.
%
% This routine should be checked.
%
% INPUT
%      strike: Fault strike measured CW from N (deg.)
%         dip: Fault dip measured CW from horizontal (deg.)
%        rake: Slip direction measured CCW in fault plane (deg.)
%     ray_phi: Ray take-off azimuth measured CCW from E (MATLAB convention; deg.)
%   ray_theta: Ray take-off dip measured CCW from E,N-plane (MATLAB convention; deg.)
%
% OUTPUT
%          up: Relative P-wave amplitude (positive = compression)
%         ush: Relative Sh-wave amplitude
%         usv: Relative Sv-wave amplitude


% Get ray coordinates in fault-oriented coordinate system
[xray,yray,zray] = sph2cart(ray_phi*pi/180,ray_theta*pi/180,1);
% Rotate ray vector into fault-plane coordinates
Nray = [xray(:)'; yray(:)'; zray(:)'];
Nray = rotz(-rake)*rotx(dip-180)*rotz(strike-90)*Nray; % The 180 shift to the dip is because equations are for hanging wall
% Ray polar coordinates in fault-aligned coordinate system
[phi,theta] = cart2sph(Nray(1,:),Nray(2,:),Nray(3,:));

% Compute double couple displacements
up  = sin(2*theta).*cos(phi); % P-wave
ush = cos(theta).*sin(phi); % Horizontally polarized S-wave
usv = cos(2*theta).*cos(phi); % Vertically polarized S-wave
