function [rhseta,rhsu] = SWERHS1D(eta,u,time)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

% form flux differences at faces
df = zeros(Nfp*Nfaces,K);
df(:) = 0.5*g*nx(:).*(eta(vmapM)-eta(vmapP));


dg = zeros(Nfp*Nfaces,K);
dg(:) = 0.5*h0*nx(:).*(u(vmapM)-u(vmapP));
% compute right hand sides of the semi-discrete PDE


rhseta = -h0*rx.*(Dr*u) + LIFT*(Fscale.*(dg));
rhsu = -g*rx.*(Dr*eta) + LIFT*(Fscale.*(df));
return
