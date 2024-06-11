function [u] = Burgers1DSSP(u,epsilon,xL,xR,FinalTime)

% function [u] = Burgers1D(u,nu,xL,xR,FinalTime)
% Purpose  : Integrate 1D Burgers equation until 
%            FinalTime starting with
%            initial condition u in the domain [xL,xR].
%            sqrt(epsilon) is the coefficient of viscosity

Globals1D;
time = 0;

% Runge-Kutta residual storage  
resu = zeros(Np,K); 

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
CFL=0.25; umax = max(max(abs(u)));
dt = CFL* min(xmin/umax,xmin^2/sqrt(epsilon));
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps; 

F = Filter1D(N,2,16);
% outer time step loop 
for tstep=1:Nsteps
    v1=u+dt*BurgersRHS1DExact(u,epsilon,xL,xR,time);
    v2=0.25*(3*u+v1+dt*BurgersRHS1DExact(v1,epsilon,xL,xR,time+dt));
    u=1/3*(u+2*v2+2*dt*BurgersRHS1DExact(v2,epsilon,xL,xR,time+0.5*dt));
    % Increment time
    time = time+dt
    u=F*u;
end;

time

return
