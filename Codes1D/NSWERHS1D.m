function [rhsh,rhsq] = NSWERHS1D(h,q,time,bx)
% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;
u=q./h;

% Define field differences at faces
dh = zeros(Nfp*Nfaces,K); dh(:) = h(vmapM)-h(vmapP);
dq = zeros(Nfp*Nfaces,K); dq(:) = q(vmapM)-q(vmapP);
% Evaluate nonlinear flux
dfq = zeros(Nfp*Nfaces,K); dfq(:) = q(vmapM).*u(vmapM)+0.5*g*h(vmapM).^2-(q(vmapP).*u(vmapP)+0.5*g*h(vmapP).^2);


%Impose boundary conditions
qin=2.5*1; uin=2.5;  hin=1;
qout=0*0.1; uout=0;  hout=0.1;




dh(mapI)=(h(vmapI)-hin);
dh(mapO)=(h(vmapO)-hout);

dq(mapI)=(q(vmapI)-qin);
dq(mapO)=(q(vmapO)-qout);

dfq(mapI) = q(vmapI)*u(vmapI)+0.5*g*h(vmapI).^2-(qin*uin+0.5*g*hin.^2);
dfq(mapO) = q(vmapO)*u(vmapO)+0.5*g*h(vmapO).^2-(qout*uout+0.5*g*hout.^2);


%LF
maxvel = max(max(abs([u+sqrt(g*h),u-sqrt(g*h)])));

% flux term
df = nx.*(dq/2.0 ) - maxvel/2.0.*dh;
dg = nx.*(dfq/2.0 ) - maxvel/2.0.*dq;


% compute right hand sides of the semi-discrete PDE
rhsh = -rx.*(Dr*q) + LIFT*(Fscale.*(df));
rhsq = -rx.*(Dr*(q.*u+0.5*g*h.^2)) + LIFT*(Fscale.*(dg))-g*h.*bx;
return
