close all
% Driver script for solving the 1D advection equations
Globals1D;
A=1;
l=2*pi;
T=2*pi;
omega=2*pi/T;
k=2*pi/l;
g=9.81;
h0=omega^2/(k^2*g);

% Order of polymomials used for approximation 
N = 6;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(0,2*pi,20);

% Initialize solver and construct grid and metric
StartUp1D;
% vmapP(1) = vmapM(end);
% vmapP(end) = vmapM(1);

% Set initial conditions
eta = A*cos(omega*0-k*x);
u= omega/(k*h0)*eta;

FinalTime = 5;
[eta,u] = SWE1D(eta,u,FinalTime);
plot(x(:),eta(:))
hold on


etaa= A*cos(omega*FinalTime-k*x);
ua=omega/(k*h0)*etaa;
plot(x(:),etaa(:))

figure
plot(x(:),u(:))
hold on
plot(x(:),ua(:))
size(x)

err = ua - u; % compute point-wise error
M = inv(V*V'); % mass matrix
errL2 = zeros(K,1);
for k = 1 : K
    errL2(k) = err(:,k)'*diag(J(:,k))*M*err(:,k);
end
errL2 = sqrt(sum(errL2)) % Global L^2-norm of error


