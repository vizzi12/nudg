function [rhsu] = RotatingHillLShapeAdvecRHS2DUpwind(u,time,cx,cy, G)
% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

C = max(abs(G.nx(:).*cx(G.vmapM) + G.ny(:).*cy(G.vmapM)));
% form field differences at faces
du = zeros(G.Nfp*G.Nfaces, G.K); 
du(:) = (G.nx(:).*cx(G.vmapM) + G.ny(:).*cy(G.vmapM)).*(u(G.vmapM) - u(G.vmapP)) - (0.5*(G.nx(:).*cx(G.vmapM) + G.ny(:).*cy(G.vmapM)).*(u(G.vmapM) - u(G.vmapP)) + 0.5*C*(u(G.vmapM) - u(G.vmapP)));

% impose boundary condition
lam = 1/8;
uin = Advec2Dsol(G.x(G.vmapB),G.y(G.vmapB),time,lam);
du(G.mapB) = (G.nx(G.mapB).*cx(G.vmapB) + G.ny(G.mapB).*cy(G.vmapB)).*(u(G.vmapB) - uin) - (0.5*(G.nx(G.mapB).*cx(G.vmapB) + G.ny(G.mapB).*cy(G.vmapB)).*(u(G.vmapB) - uin) + 0.5*C*(u(G.vmapB) - uin));

% compute right hand sides of the semi-discrete PDE
[ux, uy] = Grad2D(u, G);

rhsu = -cx.*ux - cy.*uy + G.LIFT*(G.Fscale.*(du));
return