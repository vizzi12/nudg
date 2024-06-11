function [U] =  analyticBC2D(xin, yin, nxin, nyin, mapW, mapD,mapO,U,t)
  
%  function [Q] = ForwardStepBC2D(xin, yin, nxin, nyin, mapI, mapO, mapW, mapC, Q, time);
% Purpose: Impose channel boundary conditions on 2D Euler equations on weak form

% extract conserved variables
h = U(:,:,1); q = U(:,:,2); p = U(:,:,3);

% Inflow conditions -- uniform inflow


h(mapD) = cos(-t + xin(mapD)).*cos(-t + yin(mapD)) + 2;
q(mapD) = sin(xin(mapD)-t);
p(mapD) = sin(yin(mapD)-t);


U(:,:,1)=h; U(:,:,2)=q; U(:,:,3)=p;
return