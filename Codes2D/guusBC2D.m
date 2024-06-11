function [U] =  guusBC2D(xin, yin, nxin, nyin, mapW, mapD,mapO,U,t)
  
%  function [Q] = ForwardStepBC2D(xin, yin, nxin, nyin, mapI, mapO, mapW, mapC, Q, time);
% Purpose: Impose channel boundary conditions on 2D Euler equations on weak form
g=9.81;
% extract conserved variables
h = U(:,:,1); q = U(:,:,2); p = U(:,:,3);

% Inflow conditions -- uniform inflow

idx=(xin(mapD) < 0);
h(mapD(idx))= -h(mapD(idx))+2*1;
q(mapD(idx)) = -q(mapD(idx))+2*(3*sqrt(g));
p(mapD(idx))=-p(mapD(idx));

% 
% qO = q(mapO); pO = p(mapO); 
% nxO = nxin(mapO);   nyO = nyin(mapO);
% %outflow
% q(mapO) = -qO + 2*nxO.*(nxO.*qO + nyO.*pO);
% p(mapO) = -pO + 2*nyO.*(nxO.*qO + nyO.*pO);

% q(mapO) =  nxO.^2.*qO+2*nxO.*nyO.*pO-nyO.^2.*qO;
% p(mapO) = -nxO.^2.*pO+2*nxO.*nyO.*qO+nyO.^2.*pO;
% h(mapO)=-h(mapO)+2.0;

% q(mapO)=2*nxO.^2.*qO+nxO.*nyO.*pO+nyO.^2.*qO;
% p(mapO)=nxO.^2.*pO+nxO.*nyO.*qO+2*nyO.^2.*pO;

% q(mapO)=nxO.*(nxO.*qO-nyO.*pO);
% p(mapO)=nyO.*(nxO.*qO-nyO.*pO);

% q(mapO)=0;
% p(mapO)=0;

% Wall conditions -- reflective, isothermal, i.e., n.u=0, T=T(t=0)
qW = q(mapW); pW = p(mapW); 
nxW = nxin(mapW);   nyW = nyin(mapW);

% reverse flow in normal direction in ghost elements
q(mapW) = qW - 2*nxW.*(nxW.*qW + nyW.*pW);
p(mapW) = pW - 2*nyW.*(nxW.*qW + nyW.*pW);

U(:,:,1)=h; U(:,:,2)=q; U(:,:,3)=p;
return