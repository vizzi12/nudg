% Driver script for solving the 1D advection equations with variable coefficient
clear all
close all
clc


Globals1D;

% Polynomial order used for approximation 
N = 4;

% Read in Mesh
[Nv, VX, K, EToV] = MeshGen1D(-1,1,20);

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
u = sin(x*pi);

% Solve Problem
FinalTime = 0.8;
[u,time] = Heat1D(u,FinalTime);

plot(x,u)
hold on
plot(x,exp(-pi^2*FinalTime)*sin(pi*x),'*')
