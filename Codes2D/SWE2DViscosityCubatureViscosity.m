function [U] = SWE2DViscosityCubatureViscosity(U,FinalTime,BC,BC2,S,fluxtype,viscositytype)
% function [u,time] = Advec2D(u, FinalTime, cx, cy)
% Purpose : Integrate 2D advection equation until FinalTime starting with
% initial condition u

time = 0;
% build cubature information
CubatureOrder = floor(2*(S.N+1)*3/2); cub = CubatureVolumeMesh2D(CubatureOrder,S);
% build Gauss node data
NGauss = floor((S.N+1)*2); gauss = GaussFaceMesh2D(NGauss,S);

S.dtscale = dtscale2D(S);
% compute time step size
rLGL = JacobiGQ(0,0,S.N); rmin = abs(rLGL(1)-rLGL(2));
dt = min(S.dtscale)*rmin*2/3; dt=dt/80;
tstep = 0;

S2 = viscositySmoothStartUp(S);

% hsave=zeros([size(U(:,:,1)),ceil(FinalTime/10/dt)]);
% count=1;

while (time<FinalTime)
    tstep= tstep+1;
    if(time+dt>FinalTime), dt = FinalTime-time; end
    rhsU=SWEViscosityRHS2DCubatureViscosity(U,BC,BC2,fluxtype,time,S,S2,gauss,cub,viscositytype);
    U1 = U + dt*rhsU;
          
    rhsU=SWEViscosityRHS2DCubatureViscosity(U1,BC,BC2,fluxtype,time,S,S2,gauss,cub,viscositytype);
    U = (U + U1 + dt*rhsU)/2;

    % Increment time
    time = time+dt
%     dt = SWEDT2D(U,S);
    if mod(tstep,20)==0
        PlotField2D(S.N, S.x, S.y, U(:,:,1),S,1);
%         PlotField2D(S.N, S.x, S.y, U(:,:,1)+(8<= S.x & S.x<=12).*(0.2-0.05*(S.x-10).^2),S,1);
%         zlim([0 1])
        drawnow
%         hsave(:,:,count)=U(:,:,1);
%         count=count+1;
    end
end

return