function [rhsu,rhsv] = WaveRHS1D(u,v,time)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

% form field differences at faces
du = zeros(Nfp*Nfaces,K); 
du(:) = (u(vmapM)-u(vmapP)).*nx(:)*0.5;

dv = zeros(Nfp*Nfaces,K); 
dv(:) = (v(vmapM)-v(vmapP)).*nx(:)*0.5;

% impose boundary condition at x=0
uin = cos(-pi+time);
uout = cos(pi+time);
du (mapI) = (u(vmapI)-uin).*nx(mapI)*0.5;
du (mapO) = (u(vmapO)-uout).*nx(mapO)*0.5;
% du (mapO) = 0;
% du (mapI) = 0;

vin = sin(-pi+time);
vout = sin(pi+time);
dv (mapI) = (v(vmapI)-vin).*nx(mapI)*0.5;
% dv (mapI) = 0;
dv (mapO) = (v(vmapO)-vout).*nx(mapO)*0.5;
% dv (mapO) = 0;

% compute right hand sides of the semi-discrete PDE
rhsu = -rx.*(Dr*v) + LIFT*(Fscale.*(dv));%-sin(-x + time) + cos(-x + time);
rhsv = -rx.*(Dr*u) + LIFT*(Fscale.*(du));%-cos(-x + time) + sin(-x + time);
return
