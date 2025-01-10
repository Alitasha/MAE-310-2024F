clear;
close all;
clc;

%% Note

%% Parameters
traction = 10000.0;       % Pa or N/m^2
RR        = 0.5;                % m
LL         = 4.0;               % m
radius   = 1;                 % m, variable
theta     = 0;                 % rad, variable
E           = 1e9;              % Pa or N/m^2
Poisson = 0.3;               % unity
u = E/(1+2*Poisson);     % Pa or N/m^2, shear Modulus
lamda = ( Poisson * E ) / ( (1+Poisson) * (1-2*Poisson) ); % 拉梅系数
D = [lamda+2*u lamda 0; lamda lamda+2*u 0; 0 0 u];  % Modulus Matrix

%% Exact solution for stress
sigma_rr = @(r, theta)   (traction/2) * (1 - RR^2 / radius^2) + (traction / 2) * (1 - 4 * RR^2 / radius^2 + 3 * RR^4 / radius^4) * cos(2 * theta);
sigma_tt = @(r, theta)   (traction/2) * (1 + RR^2 / radius^2) - (traction / 2) * (1 + 3 * RR^4 / radius^4) * cos(2 * theta);
sigma_rt = @(r, theta)   -(traction/2) * (1 + 2 * R^2 / radius^2 - 3 * R^4 / radius^4) * sin(2 * theta);

%% Import mesh
PlateWithHole;
[row, col] = size(msh.QUADS);
IEN= msh.QUADS(:,1:col-1);  %% 维数：单元数 * 4 （2048*4）
nel = row;

