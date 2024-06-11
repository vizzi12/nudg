function [W1,W2] =  analyticViscosityBC2D(xin, yin, nxin, nyin, mapW, mapD,mapO,W1,W2,t)
  
%  function [Q] = ForwardStepBC2D(xin, yin, nxin, nyin, mapI, mapO, mapW, mapC, Q, time);
% Purpose: Impose channel boundary conditions on 2D Euler equations on weak form

% extract conserved variables
hx = W1(:,:,1); qx = W1(:,:,2); px= W1(:,:,3);
hy = W2(:,:,1); qy = W2(:,:,2); py= W2(:,:,3);

% Inflow conditions -- uniform inflow

% idx=(xin(mapD) < 12.5);
% q(mapD(idx)) = -q(mapD(idx))+2*0.18;
% h(mapD(~idx))= -h(mapD(~idx))+2*0.33; 

hx(mapD)=sin(-xin(mapD) + t).*cos(-yin(mapD) + t);
qx(mapD)=cos(-xin(mapD) + t);
px(mapD)=0;

hy(mapD)=sin(-yin(mapD) + t).*cos(-xin(mapD) + t);
qy(mapD)=0;
py(mapD)=cos(-yin(mapD) + t);


W1(:,:,1)=hx; W1(:,:,2)=qx; W1(:,:,3)=px;
W2(:,:,1)=hy; W2(:,:,2)=qy; W2(:,:,3)=py;
return