function [W1,W2] =  parabolicDerivativeBC2D(xin, yin, nxin, nyin, mapW, mapD,mapO,W1,W2,t)
  
%  function [Q] = ForwardStepBC2D(xin, yin, nxin, nyin, mapI, mapO, mapW, mapC, Q, time);
% Purpose: Impose channel boundary conditions on 2D Euler equations on weak form

% extract conserved variables
hx = W1(:,:,1); qx = W1(:,:,2); px= W1(:,:,3);
hy = W2(:,:,1); qy = W2(:,:,2); py= W2(:,:,3);

% Inflow conditions -- uniform inflow

% idx=(xin(mapD) < 12.5);
% q(mapD(idx)) = -q(mapD(idx))+2*0.18;
% h(mapD(~idx))= -h(mapD(~idx))+2*0.33;


% hy(mapD)=-hy(mapD);
% hy(mapW)=-hy(mapW);
% 
% qy(mapD)=-qy(mapD);
% qy(mapW)=-qy(mapW);
% 
% py(mapD)=-py(mapD);
% py(mapW)=-py(mapW);
% 
% px(mapD)=-px(mapD);
% px(mapW)=-px(mapW);
% 
% qx(mapD)=-qx(mapD);
% qx(mapW)=-qx(mapW);
% 
% hx(mapD)=-hx(mapD);
% hx(mapW)=-hx(mapW);


% Wall conditions -- reflective, isothermal, i.e., n.u=0, T=T(t=0)
hxW = hx(mapW); hyW=hy(mapW);
qxW = qx(mapW); qyW=qy(mapW); 
pxW = px(mapW); pyW=py(mapW);
 
nxW = nxin(mapW);   nyW = nyin(mapW);
 
hx(mapW) = hxW - 2*nxW.*(nxW.*hxW + nyW.*hyW);
hy(mapW) = hyW - 2*nyW.*(nxW.*hxW + nyW.*hyW);

% reverse flow in normal direction in ghost elements
qx(mapW) = nxW.^4.*qxW + 2*nxW.^3.*nyW.*pxW + 2*nxW.^2.*nyW.^2.*pyW - 2*nxW.*nyW.^3.*qyW + nyW.^4.*qxW;
px(mapW) = -nxW.^4.*pxW + 2*nyW.*(qxW-pyW).*nxW.^3 + 2*nyW.^2.*(qyW+pxW).*nxW.^2 + nyW.^4.*pxW;

qy(mapW) = nxW.^4.*qyW + 2*nyW.^2.*(qyW+pxW).*nxW.^2 - 2*nyW.^3.*(qxW-pyW).*nxW - nyW.^4.*qyW; 
py(mapW) = nxW.^4.*pyW - 2*nxW.^3.*nyW.*pxW + 2*nxW.^2.*nyW.^2.*qxW + 2*nxW.*nyW.^3.*qyW + nyW.^4.*pyW;


W1(:,:,1)=hx; W1(:,:,2)=qx; W1(:,:,3)=px;
W2(:,:,1)=hy; W2(:,:,2)=qy; W2(:,:,3)=py;
return