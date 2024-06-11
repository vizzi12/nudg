function [rhseta,rhsu] = SWERHS1Dupwind(eta,u,time)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

% form flux differences at faces
deta = eta(vmapM)-eta(vmapP); du = u(vmapM)-u(vmapP);
% form flux differences at faces
df = zeros(Nfp*Nfaces,K); df(:) = 0.5*h0*nx(:).*du-0.5*sqrt(g*h0)*deta;
dg = zeros(Nfp*Nfaces,K); dg(:) = 0.5*g*nx(:).*deta-0.5*sqrt(g*h0)*du;

% Impose left outer boundary
x0 = min(x(:));
etaP=A*cos(omega*time-k*x0); uP=omega/(k*h0)*etaP;
etaM=eta(vmapI); uM=u(vmapI);
etas = 0.5*( etaP + etaM -h0/sqrt(g*h0)*(uM-uP) ); us = 0.5*( sqrt(g*h0)/h0*(etaP - etaM) + (uM+uP) );
df(mapI) = 0.5*h0*nx(mapI)*(uM-us)-0.5*sqrt(g*h0)*(etaM-etas);
dg(mapI) = 0.5*g*nx(mapI)*(etaM-etas)-0.5*sqrt(g*h0)*(uM-us);

% Impose right outer boundary
x0 = max(x(:));
etaP=A*cos(omega*time-k*x0); uP=omega/(k*h0)*etaP;
etaM=eta(vmapO); uM=u(vmapO);
etas = 0.5*( etaP + etaM +h0/sqrt(g*h0)*(uM-uP) ); us = 0.5*( sqrt(g*h0)/h0*(-etaP + etaM) + (uM+uP) );
df(mapO) = 0.5*h0*nx(mapO)*(uM-us)-0.5*sqrt(g*h0)*(etaM-etas);
dg(mapO) = 0.5*g*nx(mapO)*(etaM-etas)-0.5*sqrt(g*h0)*(uM-us);


% compute right hand sides of the semi-discrete PDE
rhseta = -h0*rx.*(Dr*u) + LIFT*(Fscale.*(df));
rhsu = -g*rx.*(Dr*eta) + LIFT*(Fscale.*(dg));
return