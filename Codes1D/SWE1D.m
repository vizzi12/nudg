function [eta,u] = SWE1D(eta,u, FinalTime)

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

Globals1D;
time = 0;

% Runge-Kutta residual storage  
reseta = zeros(Np,K); 
resu = zeros(Np,K); 

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
CFL=0.75; dt   = CFL/(2*pi)*xmin; dt = .5*dt;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps; 


% outer time step loop 
for tstep=1:Nsteps
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
        [rhseta,rhsu] = SWERHS1Dupwind(eta,u, timelocal);
        reseta = rk4a(INTRK)*reseta + dt*rhseta;
        resu= rk4a(INTRK)*resu + dt*rhsu;
        eta = eta+rk4b(INTRK)*reseta;
        u = u+rk4b(INTRK)*resu;
    end
    % Increment time
    time = time+dt;

    if mod(tstep,10)==0
        figure(1)
        plot(x,eta,'k')
        axis([min(x(:)) max(x(:)) -1 1])
        drawnow
    end
end;
return
