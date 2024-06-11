function [u] = RotatingHillLShapeAdvec2D(u, FinalTime, G)

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

time = 0.0;

% Runge-Kutta residual storage  
resu = zeros(G.Np,G.K); 

% compute time step size
dt = .125/max( ((G.N+1)^2)*.5*G.Fscale(:));

% advection speed
cx = 2*pi*G.y;
cy = -2*pi*G.x;

% outer time step loop 
while (time < FinalTime)

    if (time + dt > FinalTime)
        dt = FinalTime - time;
    end

    for INTRK = 1:5
        timelocal = time + G.rk4c(INTRK)*dt;
        [rhsu] = RotatingHillLShapeAdvecRHS2DUpwind(u, timelocal, cx, cy, G);
        resu = G.rk4a(INTRK)*resu + dt*rhsu;
        u = u+G.rk4b(INTRK)*resu;
    end;
    % Increment time
    time = time+dt
end;
return
