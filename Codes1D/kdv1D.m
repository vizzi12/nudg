function [u,time] = kdv1D(u,FinalTime,c,x0)

% function [u] = Heat1D(u,FinalTime)
% Purpose  : Integrate 1D heat equation until 
%            FinalTime starting with initial condition, u.

Globals1D;
time = 0;

% Runge-Kutta residual storage  
resu = zeros(Np, K); 

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
CFL=0.05;dt   = CFL*(xmin)^3;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;

% outer time step loop 
for tstep=1:Nsteps
  for INTRK = 1:5
    timelocal = time + rk4c(INTRK)*dt;        

    % compute right hand side of 1D advection equations
    [rhsu] = kdvRHS1D(u,timelocal,c,x0);

    % initiate and increment Runge-Kutta residuals
    resu = rk4a(INTRK)*resu + dt*rhsu;  
    
    % update fields
    u = u+rk4b(INTRK)*resu;
  end
  % Increment time
  time = time+dt

%   if mod(tstep,10)==0
%         figure(2)
%         plot(x,u,'k')
%         axis([min(x(:)) max(x(:)) 0 2.2])
%         drawnow
%   end
end
return
