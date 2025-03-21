clc; clear all
tic
%%% - Add function lab -
addpath(genpath('../stenglib-master'));
addpath(genpath('../igfe'));
load('ux.mat','-mat');
load('uy.mat','-mat');

load('fx.mat','-mat');
load('fy.mat','-mat');

% Strain & Mesh Info
L = 9.0; H = 9.0; % Geomotry Para
[nodsy,nodsx] = size(ux);
% Ele Para
nelx = nodsx - 1;
nely = nodsy - 1;
% 20x20 1% rel error
nspx = 10; % Number of strain data points per element along X
nspy = 10; % Number of strain data points per element along Y
[meshInfo,spInfo] = Mesh_Sp_Info(L,H,nelx,nely,nspx,nspy);

% Map global coordinates to local coordinates
eta = 10 ^ -5;
iter_num = 100;
[spsEl,sptEl] = Global2Local(meshInfo,spInfo,eta,iter_num);
spInfo.spsEl = spsEl; spInfo.sptEl = sptEl;

% Load
N = 0.1;
coord = meshInfo.coord;
loadPoint = find(coord(:,1)==L);
F = Load_F(meshInfo,N,loadPoint);

% Strain Data
% X direction first
ux = reshape(ux',[],1);
uy = reshape(uy',[],1);
[N,dNds,dNdt] = SP_Shape_Fun(spInfo);
[locJ,detJ] = Local_InvJ(meshInfo,dNds,dNdt);
[eps11,eps22,eps12] = Dis2Strain(meshInfo,locJ,detJ,dNds,dNdt,ux,uy);
nu = 0.3; % poisson's ratio
eps = Structured_Strain(eps11,eps22,eps12,nu);

% Calculate coefficient matrix
C = Get_C(meshInfo,spInfo,N,dNds,dNdt,locJ,eps);
% Boundary RF
[left,right,up,bot] = Edge_RF(meshInfo,spInfo,eps);

% Boundary Modulus
E_load = full(Load_Modulus(meshInfo,right,F));
Fy = Load_Fy(meshInfo,spInfo,eps,E_load);
C = Fix_Assemble(meshInfo,left,C);

% Eliminate DOF
rf_tot = -C(:,loadPoint) * E_load;
fix_DOF = sort([loadPoint*2-1;2*loadPoint]);
free_node = setdiff(1:meshInfo.nNod,loadPoint);
fb = find(coord(:,2)==L | coord(:,2)==0);
fb_DOF = sort([fb*2-1;2*fb]);
free_DOF = setdiff(1:meshInfo.ndof,cat(1,fix_DOF,fb_DOF));
C = C(free_DOF,free_node);
glo_E = (C' * C) \ (C' * rf_tot(free_DOF));
free_E = full(reshape((C' * C) \ (C' * rf_tot(free_DOF)),nodsx - 1,nodsy));
E_tot = [free_E',E_load];

% True modulus
E1 = 1.0;E2 = 8.0;
beta = (log(E2/E1)) / H;
E_load_exc = E1 * exp(beta * linspace(0,H,nodsy));
toc
