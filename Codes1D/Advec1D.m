function [u] = Advec1D(u, FinalTime,a,gfun)

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

Globals1D;
time = 0;

% Runge-Kutta residual storage  
resu = zeros(Np,K); 

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
CFL=0.375; dt   = CFL/(2*pi)*xmin; 
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps; 

% outer time step loop 
for tstep=1:Nsteps
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
        [rhsu] = AdvecRHS1D(u, timelocal, a, gfun);
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u+rk4b(INTRK)*resu;
    end
    % Increment time
    time = time+dt;

    if mod(tstep,10)==0
        figure(1)
        plot(x,u,'k')
        axis([min(x(:)) max(x(:)) -1 1])
        drawnow
    end
end
return
