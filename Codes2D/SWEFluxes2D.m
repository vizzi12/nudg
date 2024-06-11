function [F,G] = SWEFluxes2D(U)
     
% function [F,G,rho,u,v,p] = EulerFluxes2D(Q, gamma)
% Purpose: evaluate primitive variables and Euler flux functions
g=9.81;

% extract conserved variables
h = U(:,:,1); q = U(:,:,2); p = U(:,:,3);

% compute primitive variables

% compute flux functions
F = zeros(size(U)); 
F(:,:,1) = q; F(:,:,2) = q.^2./h+0.5*g*h.^2; F(:,:,3) = q.*p./h;

G = zeros(size(U));
G(:,:,1) = p; G(:,:,2) = F(:,:,3); G(:,:,3) = p.^2./h+0.5*g*h.^2;
return;