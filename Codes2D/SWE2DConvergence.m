clear all
close all

addpath 'Mesh2d v23'
addpath 'distmesh'
addpath '..\Codes1D'
addpath 'serviceRoutines'
addpath '..\RiemannSolver'

clc

% Globals2D;
g=9.81;

% Polynomial order used for approximation 

S.hgrid=0.5;
NN=2:8;
err=size(NN);

for i=1:length(NN)
S.N=NN(i);

sim='analyticCircle';
switch sim
    case {'analytic'}
        [S.p,S.EToV,S.VX,S.VY,S.K,S.Nv,S.hgrid]=analyticMesh2D(S.hgrid);
        BC  = @analyticBC2D;
        IC = @analyticIC2D;
        BC2  = @emptyDerivativeBC2D;
        FinalTime = 0.1;
        S=StartUp2D(S);
        S=analyticBCtype(S);
        S=analyticSource2D(S);
        S.straight=1:S.K;
        S.curved=[];

    case {'analyticCircle'}
        [S.p,S.EToV,S.VX,S.VY,S.K,S.Nv,S.hgrid]=analyticCircleMesh2D(S.hgrid);
        BC  = @analyticBC2D;
        IC = @analyticIC2D;
        BC2  = @emptyDerivativeBC2D;
        FinalTime = 0.1;
        S=StartUp2D(S);
        S=analyticCircleBCtype(S);
        S=analyticSource2D(S);
end

%% Solve Problem
figure(10)
U=IC(S.x,S.y);
fluxtype='HLLC';

type='Cubature';
form='weak';
viscositytype='MDA';
switch type
    case {'Cubature'}
        [U] = SWE2DCubature(U,FinalTime,BC,S,fluxtype);
    case {'CubatureLimit'}
        [U] = SWE2DCubatureLimit(U,FinalTime,BC,S,fluxtype);
    case {'Limit'}
        [U] = SWE2DLimit(U,FinalTime,BC,S,fluxtype,form);
    case {'CubatureViscosity'}
        [U] = SWE2DCubatureViscosity(U,FinalTime,BC,BC2,S,fluxtype,viscositytype);
    case {'ViscosityCubatureViscosity'}
        [U] = SWE2DViscosityCubatureViscosity(U,FinalTime,BC,BC2,S,fluxtype,viscositytype);
end

err(i)=max(abs(U(:,:,1)-(cos(-FinalTime + S.x).*cos(-FinalTime + S.y) + 2)),[],'all');
end

figure
semilogy(NN,err)





