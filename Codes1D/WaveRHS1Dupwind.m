function [rhsu,rhsv] = WaveRHS1Dupwind(u,v,time)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

% form field differences at faces
du = zeros(Nfp*Nfaces,K); 
du(:) = nx(:).*(v(vmapM)-(0.5*(v(vmapM)+v(vmapP))-0.5*nx(:).*(u(vmapP)-u(vmapM))     ));

dv = zeros(Nfp*Nfaces,K); 
dv(:) = nx(:).*(u(vmapM)-(0.5*(u(vmapM)+u(vmapP))-0.5*nx(:).*(v(vmapP)-v(vmapM))     ));


uin = cos(-pi+time)+sin(-pi-time);
uout = cos(pi+time)+sin(pi-time);
vin = -cos(-pi+time)+sin(-pi-time);
vout = -cos(pi+time)+sin(pi-time);

usl=0.5*(uin+vin-v(vmapI)+u(vmapI));
vsl=0.5*(uin+vin+v(vmapI)-u(vmapI));

usr=0.5*(v(vmapO)+u(vmapO)-vout+uout);
vsr=0.5*(vout-uout+u(vmapO)+v(vmapO));

du (mapI) =  nx(mapI).*(v(vmapI)-(0.5*(v(vmapI)+vsl)-0.5*nx(mapI).*(usl-u(vmapI))     ));
du (mapO) = nx(mapO).*(v(vmapO)-(0.5*(v(vmapO)+vsr)-0.5*nx(mapO).*(usr-u(vmapO))     ));

dv (mapI) =  nx(mapI).*(u(vmapI)-(0.5*(u(vmapI)+usl)-0.5*nx(mapI).*(vsl-v(vmapI))     ));
dv (mapO) = nx(mapO).*(u(vmapO)-(0.5*(u(vmapO)+usr)-0.5*nx(mapO).*(vsr-v(vmapO))     ));


% compute right hand sides of the semi-discrete PDE
rhsu = -rx.*(Dr*v) + LIFT*(Fscale.*(du));
rhsv = -rx.*(Dr*u) + LIFT*(Fscale.*(dv));
return
