function dt =SWEDT2D(U,S)

% function dt = EulerDT2D(Q, gamma)
% purpose: compute the time step dt for the compressible Euler equations
g=9.81;

h = U(:,:,1); q = U(:,:,2); p = U(:,:,3);
h = h(S.vmapM); q = q(S.vmapM); p = p(S.vmapM);

u = q./h; v = p./h;
c = sqrt(g*h);

dt = 1/max( ((S.N+1)^2)*.5*S.Fscale(:).*(sqrt ( u(:).^2 + v(:).^2 ) + c(:)));
return