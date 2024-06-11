% Driver script for solving the 1D advection equations with variable coefficient
clear all
close all
clc


Globals1D;

% Polynomial order used for approximation 
N = 4;

% Read in Mesh
[Nv, VX, K, EToV] = MeshGen1D(-20,30,50);

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
c=[2,1];
x0=[-10,-5];
u = c(1)/2*sech(0.5*sqrt(c(1))*(x-c(1)*0-x0(1))).^2+c(2)/2*sech(0.5*sqrt(c(2))*(x-c(2)*0-x0(2))).^2;
figure
plot(x,u)

% Solve Problem
FinalTime = 5;
[u,time] = kdv1D(u,FinalTime,c,x0);

%%
figure
plot(x,u)
hold on
ua=c(1)/2*sech(0.5*sqrt(c(1))*(x-c(1)*FinalTime-x0(1))).^2+c(2)/2*sech(0.5*sqrt(c(2))*(x-c(2)*FinalTime-x0(2))).^2;
plot(x,ua,'*')
