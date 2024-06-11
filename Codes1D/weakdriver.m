% Driver script for solving the 1D advection equations
Globals1D;

% Order of polymomials used for approximation 
N = 8;
k=160;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(0.0,160.0,k);

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
u = sin(x);

FinalTime = 26.5;
[u] = Advec1Dweak(u,FinalTime);
plot(x(:),u(:)); 

ua= sin(x-2*pi*26.5);

err = ua - u; % compute point-wise error
M = inv(V*V'); % mass matrix
errL2 = zeros(K,1);
for k = 1 : K
    errL2(k) = err(:,k)'*diag(J(:,k))*M*err(:,k);
end
errL2 = sqrt(sum(errL2)) % Global L^2-norm of error



