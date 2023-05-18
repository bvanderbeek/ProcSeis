%% Test AxiSEM 1D Model
clear;
close all;
clc;

% Input
% Grid vectors (km)
x = (-1500:20:1500)';
z = -(0:20:1000)';

% Create 2D mesh
[Zm,Xm] = meshgrid(z,x);
[nx,nz] = size(Xm);

% Read in values from AxiSEM 1D model
[rho,vp,vs] = sf_get_axisem_1D(Xm(:),zeros(nx*nz,1),Zm(:),'');

% Reshape
rho = reshape(rho,nx,nz);
vp  = reshape(vp,nx,nz);
vs  = reshape(vs,nx,nz);

figure;
imagesc(x,z,vs');
set(gca,'ydir','normal');
axis image;
colormap(jet);
