function [U] = SWE2DLimit(U,FinalTime,BC,S,fluxtype,form)
% function [u,time] = Advec2D(u, FinalTime, cx, cy)
% Purpose : Integrate 2D advection equation until FinalTime starting with
% initial condition u
% Globals2D;
time = 0;

% compute time step size
rLGL = JacobiGQ(0,0,S.N); rmin = abs(rLGL(1)-rLGL(2));
dtscale = dtscale2D(S); dt = min(dtscale)*rmin*2/3; dt=dt/20;
tstep = 0;

% dt = 1.5*SWEDT2D(U,S);

limit = limitStart(S);
[U] = SWELimiter2Dn(U,BC,S,limit,0);

while (time<FinalTime)
    tstep= tstep+1;
    if(time+dt>FinalTime), dt = FinalTime-time; end
        rhsU=SWERHS2D(U,BC,fluxtype,form,time,S);
        U1 = U + dt*rhsU;
        U1 = SWELimiter2Dn(U1,BC,S,limit,time+dt);
          
        rhsU=SWERHS2D(U1,BC,fluxtype,form,time+dt,S);
        U = (U + U1 + dt*rhsU)/2;
        U = SWELimiter2Dn(U,BC,S,limit,time+dt);

        % Increment time
        time = time+dt
%         dt = SWEDT2D(U,S);
%     if mod(tstep,1)==0
%         plotSol(h,10)
%         drawnow
%     end
end

return