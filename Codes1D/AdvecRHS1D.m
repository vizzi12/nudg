function [rhsu] = AdvecRHS1D(u,time, a,gfun)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

% form field differences at faces
% alpha=0;
du = zeros(Nfp*Nfaces,K);
alpha = 0; % upwind
du(:) = 0.5*a(vmapM).*(u(vmapM)-u(vmapP)).*(nx(:)-(1-alpha));
% impose boundary condition at x=-2
uin = sin(pi*(x(1,1)-time));
du(mapI) = 0.5*a(vmapI).*(u(vmapI)-uin).*(nx(mapI)-(1-alpha));
du(mapO) = 0;rx


rhsg = gfun(x,time);

% compute right hand sides of the semi-discrete PDE
rhsu = a.*(-rx.*(Dr*u))+(LIFT*(Fscale.*du))+rhsg;

% rhsu = a.*(-rx.*(Dr*u))+(LIFT*(Fscale.*du));

return
