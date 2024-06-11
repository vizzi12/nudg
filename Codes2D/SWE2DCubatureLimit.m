function [U] = SWE2DCubatureLimit(U,FinalTime,BC,S,fluxtype)
% function [u,time] = Advec2D(u, FinalTime, cx, cy)
% Purpose : Integrate 2D advection equation until FinalTime starting with
% initial condition u
% Globals2D;
time = 0;
% build cubature information
CubatureOrder = floor(2*(S.N+1)*3/2); cub = CubatureVolumeMesh2D(CubatureOrder,S);
% build Gauss node data
NGauss = CubatureOrder+1; gauss = GaussFaceMesh2D(NGauss,S);

% compute time step size
rLGL = JacobiGQ(0,0,S.N); rmin = abs(rLGL(1)-rLGL(2));
dtscale = dtscale2D(S); dt = min(dtscale)*rmin*2/3; dt=dt/80;
tstep = 0;

% dt = SWEDT2D(U,S);

limit = limitStart(S);
[U] = SWELimiter2Dn(U,BC,S,limit,time);

while (time<FinalTime)
    tstep= tstep+1;
    if(time+dt>FinalTime), dt = FinalTime-time; end
        rhsU=SWERHS2DCubature(U,BC,fluxtype,time,S,gauss,cub);
        U1 = U + dt*rhsU;
        U1= SWELimiter2Dn(U1,BC,S,limit,time);
          
        rhsU=SWERHS2DCubature(U1,BC,fluxtype,time+dt,S,gauss,cub);
        U = (U + U1 + dt*rhsU)/2;
        U = SWELimiter2Dn(U,BC,S,limit,time+dt);

        % Increment time
        time = time+dt
%         dt = SWEDT2D(U,S);
%     if mod(tstep,100)==0
%         PlotField2D(S.N, S.x, S.y, U(:,:,1),S,1);
%         PlotField2D(S.N, S.x, S.y, U(:,:,1)+(8<= S.x & S.x<=12).*(0.2-0.05*(S.x-10).^2),S,1);
%         zlim([0 1])
%         drawnow
%     end
end

return