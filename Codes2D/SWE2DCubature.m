function [U] = SWE2DCubature(U,FinalTime,BC,S,fluxtype)
% function [u,time] = Advec2D(u, FinalTime, cx, cy)
% Purpose : Integrate 2D advection equation until FinalTime starting with
% initial condition u
% Globals2D;
time = 0;
% build cubature information
CubatureOrder = floor(2*(S.N+1)*3/2); cub = CubatureVolumeMesh2D(CubatureOrder,S);
% build Gauss node data
% NGauss = floor((S.N+1)*2);
NGauss=CubatureOrder+1;
gauss = GaussFaceMesh2D(NGauss,S);

% compute time step size
rLGL = JacobiGQ(0,0,S.N); rmin = abs(rLGL(1)-rLGL(2));
dtscale = dtscale2D(S); dt = min(dtscale)*rmin*2/3; dt=dt/80;
tstep = 0;

% dt = SWEDT2D(U,S);
% F=Filter2D(S.N,0,2,S);

% [TRI,xout,yout,interp,S] = plotStart(S.N,S);

while (time<FinalTime)
    tstep= tstep+1;
    if(time+dt>FinalTime), dt = FinalTime-time; end
        rhsU=SWERHS2DCubature(U,BC,fluxtype,time,S,gauss,cub);
        U1 = U + dt*rhsU;
%         U1=pagemtimes(F,U1);
          
        rhsU=SWERHS2DCubature(U1,BC,fluxtype,time+dt,S,gauss,cub);
        U = (U + U1 + dt*rhsU)/2;
%         U=pagemtimes(F,U);

        % Increment time
        time = time+dt
%         dt = SWEDT2D(U,S);
%     if mod(tstep,100)==0
%           PlotField2D(xout, yout, U(:,:,1),1,interp,TRI);
%         PlotField2D(S.N, S.x, S.y, U(:,:,1)+(8<= S.x & S.x<=12).*(0.2-0.05*(S.x-10).^2),S,1);
%         zlim([0 1])
%         drawnow
%     end
end

% resU = zeros(size(U));
% while (time<FinalTime)
%     tstep= tstep+1;
%     if(time+dt>FinalTime), dt = FinalTime-time; end
%     for INTRK = 1:5
%         timelocal = time + S.rk4c(INTRK)*dt;
%         rhsU=SWERHS2DCubature(U,BC,fluxtype,timelocal,S,gauss,cub);
%         resU = S.rk4a(INTRK)*resU + dt*rhsU;
%         U = U+S.rk4b(INTRK)*resU;
%     end
%     % Increment time
%     time = time+dt
% 
%     if mod(tstep,100)==0
%           PlotField2D(xout, yout, U(:,:,1),1,interp,TRI);
% %         PlotField2D(S.N, S.x, S.y, U(:,:,1)+(8<= S.x & S.x<=12).*(0.2-0.05*(S.x-10).^2),S,1);
% %         zlim([0 1])
%         drawnow
%     end
% end

return