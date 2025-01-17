% Driver script for solving the 1D advection equations
close all
clear all
clc

Globals1D;

% Order of polymomials used for approximation 
N = 4;
k=120;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(-2,2,k);

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
a=ones(size(x));
a(abs(x) <= 0.5) = 1.5;
a(N+1,45)=1;
a(1,76)=1;

% setup rhs source function g(x,t)
 chooseg = 0;
if chooseg
% syms u(x,t)
% syms a(x) 
% u(x,t) = sin(pi*(x-a(x)*t));
% ut = diff(u,t)
% ut(x, t) = -pi*cos(pi*(x - t*a(x)))*a(x)
% ux = diff(u,x)
% ux(x, t) = -pi*cos(pi*(x - t*a(x)))*(t*diff(a(x), x) - 1) 
% syms g(x,t)
% g(x,t) = ut(x,t) + a(x)*ux(x,t)
% g(x, t) = - pi*cos(pi*(x - t*a(x)))*a(x) - pi*cos(pi*(x - t*a(x)))*a(x)*(t*diff(a(x), x) - 1)

% source term (remark, equal to zero when a'(x) = 0 )
gfun = @(x,t) -pi*cos(pi*(x-a*t)).*a - a.*( pi*cos(pi*(x-a*t)).*(0 - 1)  ) ;

else
% syms u(x,t)
% syms a(x) 
% c = 1;
% u(x,t) = sin(pi*(x-c*t));
% ut = diff(u,t)
% ux = diff(u,x)
% syms g(x,t)
% g(x,t) = ut(x,t) + a(x)*ux(x,t)
% ut(x, t) = -pi*cos(pi*(t - x))
% ux(x, t) = pi*cos(pi*(t - x))
%g(x, t) = pi*cos(pi*(t - x))*a(x) - pi*cos(pi*(t - x))
gfun = @(x,t) pi*cos(pi*(t-x)).*(a-1);
end

u = sin(pi*(x - a*0));
plot(x,u); 
FinalTime = 1;
[u] = Advec1D(u,FinalTime,a, gfun);
plot(x,u); 

% exact solutions
if chooseg
    ua = sin(pi*(x - a*FinalTime)); % this one is wrong
else
    c=1;
    ua = sin(pi*(x - c*FinalTime));
end

hold on
plot(x,ua,'r'); 

err = ua - u; % compute point-wise error
M = inv(V*V'); % mass matrix
errL2 = zeros(K,1);
for k = 1 : K
    errL2(k) = err(:,k)'*diag(J(:,k))*M*err(:,k);
end
errL2 = sqrt(sum(errL2)) % Global L^2-norm of errors



