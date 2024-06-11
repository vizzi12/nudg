function [h,q,p] = SWE2D(h,q,p,FinalTime,hf,qf,pf,s1,s2,s3)
% function [u,time] = Advec2D(u, FinalTime, cx, cy)
% Purpose : Integrate 2D advection equation until FinalTime starting with
% initial condition u
Globals2D;
time = 0;
% Runge-Kutta residual storage
resh = zeros(Np,K);
resq = zeros(Np,K);
resp = zeros(Np,K);
% compute time step size
rLGL = JacobiGQ(0,0,N); rmin = abs(rLGL(1)-rLGL(2));
dtscale = dtscale2D; dt = min(dtscale)*rmin*2/3; dt=dt/4;
% outer time step loop
tstep = 0;
tri = delaunay(x,y);
Filt = CutOffFilter2D(N,0.95);
while (time<FinalTime)
    tstep= tstep+1;
    if(time+dt>FinalTime), dt = FinalTime-time; end
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
%         [rhsh,rhsq,rhsp] = SWERHS2D(h,q,p,timelocal,hf,qf,pf,s1,s2,s3);
%         [rhsh,rhsq,rhsp] = SWERHS2DHLL(h,q,p,timelocal,hf,qf,pf,s1,s2,s3);
          [rhsh,rhsq,rhsp] = SWERHS2DHLLC(h,q,p,timelocal,hf,qf,pf,s1,s2,s3);
%         [rhsh,rhsq,rhsp] = SWERHS2DHLLweak(h,q,p,timelocal,hf,qf,pf,s1,s2,s3);
%         [rhsh,rhsq,rhsp] = SWERHS2Dweak(h,q,p,timelocal,hf,qf,pf,s1,s2,s3);
%         [rhsh2,rhsq2,rhsp2,temp2]= SWERHS2DAnalytic(h,q,p, timelocal);
        rhsh = Filt*rhsh; rhsq = Filt*rhsq; rhsp = Filt*rhsp;
        resh = rk4a(INTRK)*resh + dt*rhsh;
        resq= rk4a(INTRK)*resq + dt*rhsq;
        resp= rk4a(INTRK)*resp + dt*rhsp;
        h = h+rk4b(INTRK)*resh;
        q = q+rk4b(INTRK)*resq;
        p = p+rk4b(INTRK)*resp;
    end
    % Increment time
    time = time+dt
    if mod(tstep,10)==0
        figure(10)
        trisurf(tri,x,y,h)
        drawnow
    end
end
return