function [rhsu] = BurgersRHS1DExact(u,epsilon,xL,xR,time)

% function [rhsu] = BurgersRHS1D(u,epsilon,xL, xR, time)
% Purpose  : Evaluate RHS flux in 1D viscous Burgers equation

Globals1D;

[rnew,w] = JacobiGQ(0,0,ceil(1.5*N));
Vnew = Vandermonde1D(N,rnew);
unew = Vnew*(V \ u);

fhat=((unew.^2)')*(Vnew.*w);
f=V*fhat';

% f=u.^2;



% Define field differences at faces
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapM)-u(vmapP);

% impose boundary condition at x=0
uin =1./(cosh(epsilon*(xL+5.0)-time).^2)+1.0; du(mapI) = 2.0*(u(vmapI)-uin);
uout=1./(cosh(epsilon*(xR+5.0)-time).^2)+1.0; du(mapO) = 2.0*(u(vmapO)-uout);


% Evaluate nonlinear flux
du2 = zeros(Nfp*Nfaces,K); du2(:) = (f(vmapM)-f(vmapP))/2.0;

% impose boundary condition
du2(mapI)=(f(vmapI)-uin.^2)/2.0; du2(mapO)=(f(vmapO)-uout.^2)/2.0;

% Compute flux
% maxvel = max(max(abs(u)));
maxvel=0;

% penalty scaling -- See Chapter 7.2
%tau = .25*reshape(N*N./max(2*J(vmapP),2*J(vmapM)), Nfp*Nfaces, K);
% tau=0;

% flux term
% flux = nx.*(du2/2.0 - sqrt(epsilon)*dq) - maxvel/2.0.*du - sqrt(epsilon)*tau.*du;
flux = nx.*(du2/2.0 )- maxvel/2.0.*du; %- sqrt(epsilon)*tau.*du;


% local derivatives of field
% dfdr = Dr*(u.^2/2 - sqrt(epsilon)*q); 
dfdr = Dr*(f/2); 


y=epsilon * x + 5.0 * epsilon - time;
g=-2 * sinh(y) .* ((epsilon - 1)*cosh(y) .^ 2 + epsilon) ./ cosh(y) .^ 5;

% compute right hand sides of the semi-discrete PDE
rhsu = -(rx.*dfdr - LIFT*(Fscale.*flux))+g;
return
