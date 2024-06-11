function [rhsu,rhsp] = MachRHS1Dupwind(u,p,time,M)
Globals1D;


df = zeros(Nfp*Nfaces,K); 
df(:) = nx(:).* (M(vmapM).*u(vmapM)+p(vmapM)) -  (0.5*nx(:).*(M(vmapM).*(u(vmapP)+u(vmapM))+ p(vmapP)+p(vmapM))+0.5*(M(vmapM).*(p(vmapM)-p(vmapP))+u(vmapM)-u(vmapP)  ));

dg = zeros(Nfp*Nfaces,K); 
dg(:) = nx(:).* (M(vmapM).*p(vmapM)+u(vmapM)) -  (0.5*nx(:).*(u(vmapP)+u(vmapM)+ M(vmapM).*(p(vmapP)+p(vmapM))  )+0.5*(p(vmapM)-p(vmapP)+M(vmapM).*(u(vmapM)-u(vmapP))  ));

% impose boundary condition at x=0
uin = cos(-pi-time);
uout = cos(pi-time);
pin = sin(-pi-time);
pout = sin(pi-time);

df(mapI) = nx(mapI).* (M(vmapI).*u(vmapI)+p(vmapI)) -  (0.5*nx(mapI).*(M(vmapI).*(uin+u(vmapI))+ pin+p(vmapI))+0.5*(M(vmapI).*(p(vmapI)-pin)+u(vmapI)-uin  ));
dg(mapI) = nx(mapI).* (M(vmapI).*p(vmapI)+u(vmapI)) -  (0.5*nx(mapI).*(uin+u(vmapI)+ M(vmapI).*(pin+p(vmapI))  )+0.5*(p(vmapI)-pin+M(vmapI).*(u(vmapI)-uin)  ));

df(mapO) = nx(mapO).* (M(vmapO).*u(vmapO)+p(vmapO)) -  (0.5*nx(mapO).*(M(vmapO).*(uout+u(vmapO))+ pout+p(vmapO))+0.5*(M(vmapO).*(p(vmapO)-pout)+u(vmapO)-uout  ));
dg(mapO) = nx(mapO).* (M(vmapO).*p(vmapO)+u(vmapO)) -  (0.5*nx(mapO).*(uout+u(vmapO)+ M(vmapO).*(pout+p(vmapO))  )+0.5*(p(vmapO)-pout+M(vmapO).*(u(vmapO)-uout)  ));

% compute right hand sides of the semi-discrete PDE
rhsu = -rx.*M.*(Dr*u)-rx.*(Dr*p) + LIFT*(Fscale.*(df))+(M-1).*sin(-x+time)+cos(-x+time);
rhsp = -rx.*(Dr*u)-rx.*M.*(Dr*p) + LIFT*(Fscale.*(dg))+(M-1).*cos(-x+time)+sin(-x+time);
return
