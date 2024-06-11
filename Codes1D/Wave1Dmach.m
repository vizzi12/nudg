function [u,p] = Wave1Dmach(u,p, FinalTime,M)

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

Globals1D;
time = 0;

% Runge-Kutta residual storage  
resu = zeros(Np,K); 
resp = zeros(Np,K); 

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
CFL=0.75; dt   = CFL*1/2*min(VX(2:end)-VX(1:end-1))*0.5*min(r(2:end)-r(1:end-1)); dt=dt*0.5;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps; 


% outer time step loop 
for tstep=1:Nsteps
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
        [rhsu,rhsp] = MachRHS1Dupwind(u,p, timelocal,M);
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resp= rk4a(INTRK)*resp + dt*rhsp;
        u = u+rk4b(INTRK)*resu;
        p = p+rk4b(INTRK)*resp;
    end
    % Increment time
    time = time+dt;


end;
return
