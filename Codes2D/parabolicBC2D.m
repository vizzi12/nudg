function [U] =  parabolicBC2D(xin, yin, nxin, nyin, mapW, mapD,mapO,U,t)
  
%  function [Q] = ForwardStepBC2D(xin, yin, nxin, nyin, mapI, mapO, mapW, mapC, Q, time);
% Purpose: Impose channel boundary conditions on 2D Euler equations on weak form

% extract conserved variables
h = U(:,:,1); q = U(:,:,2); p = U(:,:,3);

% Inflow conditions -- uniform inflow

idx=(xin(mapD) < 12.5);
q(mapD(idx)) = -q(mapD(idx))+2*0.18;
h(mapD(~idx))= -h(mapD(~idx))+2*0.33;

p(mapD)=-p(mapD);
% p(mapW)=-p(mapW);

% Wall conditions -- reflective, isothermal, i.e., n.u=0, T=T(t=0)
qW = q(mapW); pW = p(mapW); 
nxW = nxin(mapW);   nyW = nyin(mapW);

% reverse flow in normal direction in ghost elements
q(mapW) = qW - 2*nxW.*(nxW.*qW + nyW.*pW);
p(mapW) = pW - 2*nyW.*(nxW.*qW + nyW.*pW);

U(:,:,1)=h; U(:,:,2)=q; U(:,:,3)=p;
return