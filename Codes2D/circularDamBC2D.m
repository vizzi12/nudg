function [U] =  circularDamBC2D(xin, yin, nxin, nyin, mapW, mapD,mapO,U,t)
  
%  function [Q] = ForwardStepBC2D(xin, yin, nxin, nyin, mapI, mapO, mapW, mapC, Q, time);
% Purpose: Impose channel boundary conditions on 2D Euler equations on weak form

% extract conserved variables

% Inflow conditions -- uniform inflow

h = U(:,:,1); q = U(:,:,2); p = U(:,:,3);

h(mapD)=0.5;
q(mapD)=0;
p(mapD)=0;

U(:,:,1)=h; U(:,:,2)=q; U(:,:,3)=p;
return