function [U,hsave] = SWE2DCubatureViscosity(U,FinalTime,BC,BC2,S,fluxtype,viscositytype)
% function [u,time] = Advec2D(u, FinalTime, cx, cy)
% Purpose : Integrate 2D advection equation until FinalTime starting with
% initial condition u

time = 0;
% build cubature information
CubatureOrder = floor(2*(S.N+1)*3/2); cub = CubatureVolumeMesh2D(CubatureOrder,S);
% build Gauss node data
NGauss = CubatureOrder+1; gauss = GaussFaceMesh2D(NGauss,S);

S.dtscale = dtscale2D(S);
% compute time step size
rLGL = JacobiGQ(0,0,S.N); rmin = abs(rLGL(1)-rLGL(2));
dt = min(S.dtscale)*rmin*2/3; dt=dt/100;
tstep = 0;

S.dt=dt;
S2 = viscositySmoothStartUp(S);

hsave=0;
% hsave=zeros([size(U(:,:,1)),ceil(FinalTime/25/dt)]);
% count=1;

if strcmp(viscositytype,'EV')
    S.Eold=nan*S.x;
    S.divold=nan;
    S = computeSWEEntropy(U,S);
    S.mapM3D=cat(3,S.vmapM,S.vmapM+length(S.vmapM(:)));
    S.mapP3D=cat(3,S.vmapP,S.vmapP+length(S.vmapP(:)));
    S.omega=sum(cub.W.*(cub.V*ones(size(S.x))),'all');
end

[TRI,xout,yout,interp,S] = plotStart(S.N,S);

while (time<FinalTime)
    tstep= tstep+1;
    if(time+dt>FinalTime), dt = FinalTime-time; end

    [rhsU]=SWERHS2DCubatureViscosity(U,BC,BC2,fluxtype,time,S,S2,gauss,cub,viscositytype);
    U1 = U + dt*rhsU;
          
    [rhsU]=SWERHS2DCubatureViscosity(U1,BC,BC2,fluxtype,time+dt,S,S2,gauss,cub,viscositytype);
    U = (U + U1 + dt*rhsU)/2;

    % Increment time
    time = time+dt
    if strcmp(viscositytype,'EV')
        S.Eold=S.E;
        S.divold=S.div;
        S = computeSWEEntropy(U,S);
    end
%     dt = SWEDT2D(U,S);
    if mod(tstep,2000)==0
%         subplot(1,2,1);
        PlotField2D(xout, yout, U(:,:,1),1,interp,TRI);
%         zlim([0 1])
%         subplot(1,2,2)
%         PlotField2D(S.N, S.x, S.y, mu,S,1);
%         view(2)
        drawnow
%         hsave(:,:,count)=U(:,:,1);
%         count=count+1;
    end
end

return