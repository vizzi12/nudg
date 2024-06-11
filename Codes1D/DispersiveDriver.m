% Driver script for solving the 1D advection equations with variable coefficient
clear all
close all
clc


Globals1D;

% Polynomial order used for approximation 
N = 2;

% Read in Mesh
[Nv, VX, K, EToV] = MeshGen1D(-1,1,20);

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
u = cos(x*pi);

% Solve Problem
FinalTime = 0.05;
[u,time] = Dispersive1D(u,FinalTime);

%%
figure
plot(x,u)
hold on
plot(x,cos(pi^3*FinalTime+pi*x),'*')
