function [h,q] = NSWE1D(h,q, FinalTime,bx)

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

Globals1D;
time = 0;

% Runge-Kutta residual storage  
resh = zeros(Np,K); 
resq = zeros(Np,K); 

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
CFL=0.75; dt   = CFL/(2*pi)*xmin; dt = dt/10;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps; 

% F = Filter1D(N,0,8);
% outer time step loop 
for tstep=1:Nsteps
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
%         [rhsh,rhsq] = NSWERHS1D(h,q, timelocal,bx);
%         [rhsh,rhsq] = NSWERHS1DHLL(h,q, timelocal,bx);
        [rhsh,rhsq]=NSWERHS1DHLLViscosity(h,q,time,bx);
%         [rhsh,rhsq]=NSWERHS1DHLLViscosity2(h,q,time,bx);
%         [rhsh,rhsq]=NSWERHS1DRPViscosity(h,q,time,bx);
        resh = rk4a(INTRK)*resh + dt*rhsh;
        resq= rk4a(INTRK)*resq + dt*rhsq;
        h = h+rk4b(INTRK)*resh;
        q = q+rk4b(INTRK)*resq;
    end
%     h=F*h;
%     q=F*q;
%     h = SlopeLimitN(h);
%     q = SlopeLimitN(q);


    % Increment time
    time = time+dt

    if mod(tstep,10)==0
        figure(1)
        plot(x,h,'k')
        axis([min(x(:)) max(x(:)) 0 4])
        drawnow
    end
end;
return
