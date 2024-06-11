function [rhsh,rhsq] = NSWERHS1DHLLViscosity2(h,q,time,bx)
% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;
u=q./h;

a=sqrt(g*h);

% hS=0.5*(h(vmapL)+h(vmapR))-0.25*(u(vmapR)-u(vmapL)).*(h(vmapL)+h(vmapR))./(a(vmapL)+a(vmapR));

hS=1/g*(0.5*(a(vmapL)+a(vmapR)) + 0.25 *(u(vmapL)-u(vmapR))).^2;

qL=sqrt( 0.5*(hS+h(vmapL)).*hS./h(vmapL).^2 );
qL(hS <= h(vmapL)) = 1;

qR=sqrt( 0.5*(hS+h(vmapR)).*hS./h(vmapR).^2 );
qR(hS <= h(vmapR)) = 1;

sL=u(vmapL)-a(vmapL).*qL;
sR=u(vmapR)+a(vmapR).*qR;

qin=0; uin=0;  hin=3.5;
qout=0; uout=0;  hout=1.25;

FLh= q(vmapL);
FRh= q(vmapR);

FLh(mapI)= qin;
FRh(mapO)= qout;

fh=( sR.*FLh-sL.*FRh+sR.*sL.*(h(vmapR)-h(vmapL)) ) ./ (sR - sL);
fh(sL >= 0) = FLh(sL >= 0);
fh(sR <= 0) = FRh(sR <= 0);

FLq= q(vmapL).*u(vmapL)+0.5*g*h(vmapL).^2;
FRq= q(vmapR).*u(vmapR)+0.5*g*h(vmapR).^2;

FLq(mapI)= qin.*uin+0.5*g*hin.^2;
FRq(mapO)= qout.*uout+0.5*g*hout.^2;

fq=( sR.*FLq-sL.*FRq+sR.*sL.*(q(vmapR)-q(vmapL)) ) ./ (sR - sL);
fq(sL >= 0) = FLq(sL >= 0);
fq(sR <= 0) = FRq(sR <= 0);

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

hhat=invV*h;
hhat=hhat(2:Np,:);

m=(1:N)';
b=m.^(-N)./sqrt(sum(m.^(-2*N)));
M = inv(V*V');

hL2 = zeros(K,1);
for k = 1 : K
    hL2(k) = h(:,k)'*diag(J(:,k))*M*h(:,k);
end

hhat=sqrt(  hhat.^2+ b.^2*hL2');

for i=1:N
    hhat(i,:)=max(hhat(min(i,N-1):N,:));
end

A=[ones(N,1),-log(m)];

tau=A\log(hhat);
tau=tau(2,:);

muMax=cMax*(hgrid/N)*max(abs([u+a;u-a]));

mu=tau*0+(tau<1).*muMax+(1 <= tau & tau <= 3).*(1-(tau-1)/2);

rnew = JacobiGL(0,0,2);
Vnew=Vandermonde1D(2,rnew);

mu2=zeros(3,K);

muavg=(mu(1:end-1)+mu(2:end))/2;
mu2(1,:)=[mu(1),muavg];
mu2(2,:)=mu;
mu2(3,:)=[muavg mu(K)];

muModal=Vnew\mu2;

mu=V*[muModal; zeros(N-2,K)];

gh=mu.*wh;
gq=mu.*wq;

dgh = zeros(Nfp*Nfaces,K); dgh(:)=(gh(vmapM)-gh(vmapP))/2;
dqh = zeros(Nfp*Nfaces,K); dqh(:)=(gq(vmapM)-gq(vmapP))/2;

% compute right hand sides of the semi-discrete PDE
rhsh = -rx.*(Dr*q)+rx.*(Dr*gh) + LIFT*(Fscale.*(dfh))-LIFT*(Fscale.*(nx.*dgh));
rhsq = -rx.*(Dr*(q.*u+0.5*g*h.^2))+rx.*(Dr*gq) + LIFT*(Fscale.*(dfq))-LIFT*(Fscale.*(nx.*dqh)); %-g*h.*bx
return
