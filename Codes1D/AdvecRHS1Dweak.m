function [rhsu] = AdvecRHS1Dweak(u,time, a)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

% flux
du = zeros(Nfp*Nfaces,K); 
du(:) = (u(vmapM)+u(vmapP)).*nx(:)*a*0.5;

% impose boundary condition at x=0
uin = -sin(a*time);
du (mapI) = (u(vmapI)+uin).*nx(mapI)*a*0.5;
du (mapO) = (u(vmapO)*2).*nx(mapO)*a*0.5;

% compute right hand sides of the semi-discrete PDE
rhsu = +a*rx.*((V*V')*((V*V')\Dr)'*u) - LIFT*(Fscale.*(du));
return
