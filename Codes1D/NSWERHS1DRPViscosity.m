function [rhsh,rhsq] = NSWERHS1DRPViscosity(h,q,time,bx)
% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;
u=q./h;

a=sqrt(g*h);

fh=zeros(size(vmapL));
fq=zeros(size(vmapL));
hL=h(vmapL);
hR=h(vmapR);

uL=u(vmapL);
uR=u(vmapR);

hL(mapI)=3.5;
hR(mapO)=1.25;

uL(mapI)=0;
uR(mapO)=0;

for i=1:length(vmapL)
    [hLOCAL,uLOCAL]=SWRPExact(hL(i),uL(i),0,hR(i),uR(i),0,1,0,0,1);
    fh(i)=hLOCAL*uLOCAL;
    fq(i)=hLOCAL*uLOCAL^2+0.5*g*hLOCAL^2;
end

% flux term
dfh = zeros(Nfp*Nfaces,K);
dfq = zeros(Nfp*Nfaces,K);

dfh(:) = nx(:).*(q(vmapM)-fh);
dfq(:) = nx(:).*(q(vmapM).*u(vmapM)+0.5*g*h(vmapM).^2-fq);

dq = zeros(Nfp*Nfaces,K); dq(:) = (q(vmapM)-q(vmapP))/2;
dh = zeros(Nfp*Nfaces,K); dh(:) = (h(vmapM)-h(vmapP))/2;

wh = rx.*(Dr*h)- LIFT*(Fscale.*(nx.*dh));
wq = rx.*(Dr*q)- LIFT*(Fscale.*(nx.*dq));

cBeta=0.5;
cMax=1;

muBeta=cBeta*abs(rx.*(Dr*u))*(hgrid/N)^2;
muMax=cMax*(hgrid/N)*max(abs([u+a;u-a]));

mu=min(muBeta,muMax.*ones(Np,1));

gh=mu.*wh;
gq=mu.*wq;

dgh = zeros(Nfp*Nfaces,K); dgh(:)=(gh(vmapM)-gh(vmapP))/2;
dqh = zeros(Nfp*Nfaces,K); dqh(:)=(gq(vmapM)-gq(vmapP))/2;

% compute right hand sides of the semi-discrete PDE
rhsh = -rx.*(Dr*q)+rx.*(Dr*gh) + LIFT*(Fscale.*(dfh))-LIFT*(Fscale.*(nx.*dgh));
rhsq = -rx.*(Dr*(q.*u+0.5*g*h.^2))+rx.*(Dr*gq) + LIFT*(Fscale.*(dfq))-LIFT*(Fscale.*(nx.*dqh))-g*h.*bx;
return
