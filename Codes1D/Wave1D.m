function [u,v] = Wave1D(u,v, FinalTime)

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

Globals1D;
time = 0;

% Runge-Kutta residual storage  
resu = zeros(Np,K); 
resv = zeros(Np,K); 

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
CFL=0.75; dt   = CFL/(2*pi)*xmin; dt = .5*dt;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps; 


% outer time step loop 
for tstep=1:Nsteps
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
        [rhsu,rhsv] = WaveRHS1D(u,v, timelocal);
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv= rk4a(INTRK)*resv + dt*rhsv;
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
    end
    % Increment time
    time = time+dt;


end;
return
