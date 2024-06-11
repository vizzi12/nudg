close all
% Driver script for solving the 1D advection equations
Globals1D;


% Order of polymomials used for approximation 
N = 6;
k=20;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(-pi,pi,k);

% Initialize solver and construct grid and metric
StartUp1D;
% vmapP(1) = vmapM(end);
% vmapP(end) = vmapM(1);

% Set initial conditions
u = cos(x);
v= sin(x);

FinalTime = 1.5;
[u,v] = Wave1D(u,v,FinalTime);
plot(x(:),u(:)); 
hold on


ua= cos(x-FinalTime);
plot(x(:),ua(:))

err = ua - u; % compute point-wise error
M = inv(V*V'); % mass matrix
errL2 = zeros(K,1);
for k = 1 : K
    errL2(k) = err(:,k)'*diag(J(:,k))*M*err(:,k);
end
errL2 = sqrt(sum(errL2)) % Global L^2-norm of error


