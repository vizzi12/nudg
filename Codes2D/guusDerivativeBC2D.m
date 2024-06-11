function [W1,W2] =  guusDerivativeBC2D(xin, yin, nxin, nyin, mapW, mapD,mapO,W1,W2,t)
  
%  function [Q] = ForwardStepBC2D(xin, yin, nxin, nyin, mapI, mapO, mapW, mapC, Q, time);
% Purpose: Impose channel boundary conditions on 2D Euler equations on weak form

% extract conserved variables
hx = W1(:,:,1); qx = W1(:,:,2); px= W1(:,:,3);
hy = W2(:,:,1); qy = W2(:,:,2); py= W2(:,:,3);

%outflow conditions 
% hxO = hx(mapO); hyO=hy(mapO);
% uxO = qx(mapO); uyO=qy(mapO); 
% vxO = px(mapO); vyO=py(mapO);
% nxO = nxin(mapO);   nyO = nyin(mapO);
% 
% hx(mapO) = hxO - 2*nxO.*(nxO.*hxO + nyO.*hyO);
% hy(mapO) = hyO - 2*nyO.*(nxO.*hxO + nyO.*hyO);

% qx(mapO) =  nxO .^ 4 .* uxO + 2 .* nyO .* (uyO + vxO) .* nxO .^ 3 - 2 .* nyO .^ 2 .* (uxO - 2 .* vyO) .* nxO .^ 2 - 2 .* nyO .^ 3 .* (uyO + vxO) .* nxO + nyO .^ 4 .* uxO;
% 
% px(mapO) =  -(nxO .^ 4) .* vxO + (2 .* nyO .* (uxO - vyO) .* nxO .^ 3) + 4 .* (nyO .^ 2) .* (uyO + vxO / 2) .* (nxO .^ 2) - (2 .* nyO .^ 3 .* (uxO - vyO) .* nxO) - (nyO .^ 4) .* vxO;
% 
% qy(mapO) = -nxO .^ 4 .* uyO + 2 .* nyO .* (uxO - vyO) .* nxO .^ 3 + 2 .* nyO .^ 2 .* (uyO + 2 .* vxO) .* nxO .^ 2 - 2 .* nyO .^ 3 .* (uxO - vyO) .* nxO - nyO .^ 4 .* uyO;
% 
% py(mapO) =  (nxO .^ 4) .* vyO - (2 .* nyO .* (uyO + vxO) .* nxO .^ 3) + 4 .* (nyO .^ 2) .* (uxO - vyO / 2) .* (nxO .^ 2) + (2 .* nyO .^ 3 .* (uyO + vxO) .* nxO) + (nyO .^ 4) .* vyO;

% qx(mapO) =  nxO.^4.*uxO+2*nxO.^3.*nyO.*uyO+2*nxO.^2.*nyO.^2.*vyO-2*nxO.*nyO.^3.*vxO+nyO.^4.*uxO;
% px(mapO) =  nxO.^4.*vxO+2.*nyO.^2.*(uyO+vxO).*nxO.^2-2*nyO.^3.*(uxO-vyO).*nxO-nyO.^4.*vxO;

% qy(mapO) = -nxO.^4.*uyO+2*nyO.*(uxO-vyO).*nxO.^3+2*nyO.^2.*(uyO+vxO).*nxO.^2+nyO.^4.*uyO;
% py(mapO) =  nxO.^4.*vyO-2*nxO.^3.*nyO.*uyO+2*nxO.^2.*nyO.^2.*uxO+2*nxO.*nyO.^3.*vxO+nyO.^4.*vyO;
% 
% hx(mapO) = -hxO;
% hy(mapO) = -hyO;

% qx(mapO) = -uxO;
% qy(mapO) = -uyO;

% px(mapO) = -vxO;
% py(mapO) = -vyO;

% hx(mapO) = -nyO.*(nxO.*hyO-nyO.*hxO);
% hy(mapO) =  nxO.*(nxO.*hyO-nyO.*hxO);
% 
% qx(mapO) = -nyO.*(nxO.*qyO-nyO.*qxO);
% qy(mapO) =  nxO.*(nxO.*qyO-nyO.*qxO);
% 
% px(mapO) = -nyO.*(nxO.*pyO-nyO.*pxO);
% py(mapO) =  nxO.*(nxO.*pyO-nyO.*pxO);

% hx(mapO) = hxO - 2*nxO.*(nxO.*hxO + nyO.*hyO);
% hy(mapO) = hyO - 2*nyO.*(nxO.*hxO + nyO.*hyO);
% 
% qx(mapO) = qxO - 2*nxO.*(nxO.*qxO + nyO.*qyO);
% qy(mapO) = qyO - 2*nyO.*(nxO.*qxO + nyO.*qyO);
% 
% px(mapO) = pxO - 2*nxO.*(nxO.*pxO + nyO.*pyO);
% py(mapO) = pyO - 2*nyO.*(nxO.*pxO + nyO.*pyO);

% reverse flow in normal direction in ghost elements

% qx(mapO) = -nxO.^2.*qxO - 2 * nxO.*nyO.*qyO + nyO.^2.*qxO;
% px(mapO) = -nxO.^2.*pxO - 2 * nxO.*nyO.*pyO + nyO.^2.*pxO;
% 
% qy(mapO) = nxO.^2.*qyO - 2 * nxO.*nyO.*qxO - nyO.^2.*qyO;
% py(mapO) = nxO.^2.*pyO - 2 * nxO.*nyO.*pxO - nyO.^2.*pyO;

% qx(mapO) = -nxO.^4.*qxO - 2*nyO.*(qyO+pxO).*nxO.^3 + 2*nyO.^2.*(qxO-pyO).*nxO.^2 + nyO.^4.*qxO; 
% px(mapO) = nxO.^4.*pxO - 2*nxO.^3.*nyO.*qxO - 2*nxO.^2.*nyO.^2.*qyO - 2*nxO.*nyO.^3.*pyO + nyO.^4.*pxO;
% qy(mapO) = nxO.^4.*qyO - 2*nxO.^3.*nyO.*qxO - 2*nxO.^2.*nyO.^2.*pxO - 2*nxO.*nyO.^3.*pyO + nyO.^4.*qyO;
% py(mapO) = nxO.^4.*pyO - 2*nyO.^2.*(qxO-pyO).*nxO.^2 -
% 2*nyO.^3.*(qyO+pxO).*nxO - nyO.^4.*pyO;



% Wall conditions -- reflective, isothermal, i.e., n.u=0, T=T(t=0)
hxW = hx(mapW); hyW=hy(mapW);
qxW = qx(mapW); qyW=qy(mapW); 
pxW = px(mapW); pyW=py(mapW);
 
nxW = nxin(mapW);   nyW = nyin(mapW);
% 
% hx(mapW) =- hx(mapW);
% hy(mapW) = -hy(mapW);
 
hx(mapW) = hxW - 2*nxW.*(nxW.*hxW + nyW.*hyW);
hy(mapW) = hyW - 2*nyW.*(nxW.*hxW + nyW.*hyW);
% 
% % reverse flow in normal direction in ghost elements
qx(mapW) = nxW.^4.*qxW + 2*nxW.^3.*nyW.*pxW + 2*nxW.^2.*nyW.^2.*pyW - 2*nxW.*nyW.^3.*qyW + nyW.^4.*qxW;
px(mapW) = -nxW.^4.*pxW + 2*nyW.*(qxW-pyW).*nxW.^3 + 2*nyW.^2.*(qyW+pxW).*nxW.^2 + nyW.^4.*pxW;

qy(mapW) = nxW.^4.*qyW + 2*nyW.^2.*(qyW+pxW).*nxW.^2 - 2*nyW.^3.*(qxW-pyW).*nxW - nyW.^4.*qyW; 
py(mapW) = nxW.^4.*pyW - 2*nxW.^3.*nyW.*pxW + 2*nxW.^2.*nyW.^2.*qxW + 2*nxW.*nyW.^3.*qyW + nyW.^4.*pyW;


W1(:,:,1)=hx; W1(:,:,2)=qx; W1(:,:,3)=px;
W2(:,:,1)=hy; W2(:,:,2)=qy; W2(:,:,3)=py;
return