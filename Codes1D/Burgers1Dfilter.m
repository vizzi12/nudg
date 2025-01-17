function [u] = Burgers1Dfilter(u,epsilon,xL,xR,FinalTime,s)

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


F = Filter1D(N,2,s);

% outer time step loop 
for tstep=1:Nsteps

    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
        [rhsu] = BurgersRHS1D(u,epsilon,xL,xR,timelocal);
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u+rk4b(INTRK)*resu;
    end;
    u=F*u;
    % Increment time
    time = time+dt
%     if mod(tstep,10)==0
%         figure(2)
%         plot(x,u,'k')
%         axis([min(x(:)) max(x(:)) 0 2.2])
%         drawnow
%     end
end;

time

return
