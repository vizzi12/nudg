function [rhsu] = DispersiveRHS1D(u,time)

% function [rhsu] = DispersiveLDGRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D u_xxx using LDG fluxes and periodic BC's

Globals1D;

% Define field differences at faces, incl BC
% du   = zeros(Nfp*Nfaces,K); du(:)     = u(vmapM)-u(vmapP);
% uin  = cos(pi^3*time+pi*x(vmapI));            du (mapI) = u(vmapI) - uin; 
% uout = cos(pi^3*time+pi*x(vmapO));            du (mapO) = u(vmapO) - uout;
% fluxu = nx.*du/2.0;
% 
% 
% % Compute local variable p, define differences, incl BC
% p = rx.*(Dr*u) - LIFT*(Fscale.*fluxu);
% dp   = zeros(Nfp*Nfaces,K); dp(:)    = p(vmapM)-p(vmapP);
% % pin = -pi*sin(pi^3*time + pi*x(vmapI));    dp(mapI) = p(vmapI) - pin;
% pout = -pi*sin(pi^3*time + pi*x(vmapO));    dp(mapO) = p(vmapO) - pout;
% % pin = -pi*sin(pi^3*time + pi*x(vmapI));    dp(mapI) = 0;
% % pout = -pi*sin(pi^3*time + pi*x(vmapO));    dp(mapO) = 0;
% fluxp = nx.*dp/2.0;
% 
% % Compute local variable q, define differences, incl BC
% q = rx.*(Dr*p) - LIFT*(Fscale.*fluxp);
% dq   = zeros(Nfp*Nfaces,K); dq(:)     = q(vmapM)-q(vmapP);
% qin  = -cos(pi^3*time + pi*x(vmapI))*pi^2;            dq (mapI) = q(vmapI) - qin; 
% qout = -cos(pi^3*time + pi*x(vmapO))*pi^2;            dq (mapO) = q(vmapO) - qout;
% % qin  = -cos(pi^3*time + pi*x(vmapI))*pi^2;            dq (mapI) = 0; 
% % qout = -cos(pi^3*time + pi*x(vmapO))*pi^2;            dq (mapO) = 0;
% fluxq = nx.*dq/2.0;

du   = zeros(Nfp*Nfaces,K); du(:)     = u(vmapM)-u(vmapP);
uin = u(vmapO); du (mapI) = u(vmapI) - uin;
uout = u(vmapI); du (mapO) = u(vmapO) - uout;
fluxu = nx.*du/2.0;

% Compute local variable p, define differences, incl BC
p = rx.*(Dr*u) - LIFT*(Fscale.*fluxu);
dp   = zeros(Nfp*Nfaces,K); dp(:)    = p(vmapM)-p(vmapP);
pin = p(vmapO); dp(mapI) = p(vmapI) - pin;
pout = p(vmapI); dp(mapO) = p(vmapO) - pout;
fluxp = nx.*dp/2.0;

% Compute local variable q, define differences, incl BC
q = rx.*(Dr*p) - LIFT*(Fscale.*fluxp);
dq   = zeros(Nfp*Nfaces,K); dq(:)     = q(vmapM)-q(vmapP);
qin = q(vmapO); dq (mapI) = q(vmapI) - qin;
qout = q(vmapI); dq (mapO) = q(vmapO) - qout;
fluxq = nx.*dq/2.0;

% compute right hand sides of the semi-discrete PDE
rhsu = -rx.*(Dr*q) + LIFT*(Fscale.*fluxq);
return
