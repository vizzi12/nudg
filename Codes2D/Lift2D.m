function S = Lift2D(S)

% function [LIFT] = Lift2D()
% Purpose  : Compute surface to volume lift term for DG formulation

% Globals2D;
Emat = zeros(S.Np, S.Nfaces*S.Nfp);

% face 1
faceR = S.r(S.Fmask(:,1));
V1D = Vandermonde1D(S.N, faceR); 
S.massEdge1 = inv(V1D*V1D');
Emat(S.Fmask(:,1),1:S.Nfp) = S.massEdge1;
S.V1D1=V1D;

% face 2
faceR = S.r(S.Fmask(:,2));
V1D = Vandermonde1D(S.N, faceR);
S.massEdge2 = inv(V1D*V1D');
Emat(S.Fmask(:,2),S.Nfp+1:2*S.Nfp) = S.massEdge2;
S.V1D2=V1D;

% face 3
faceS = S.s(S.Fmask(:,3));
V1D = Vandermonde1D(S.N, faceS); 
S.massEdge3 = inv(V1D*V1D');
Emat(S.Fmask(:,3),2*S.Nfp+1:3*S.Nfp) = S.massEdge3;
S.V1D3=V1D;

% inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
S.LIFT = S.V*(S.V'*Emat);
return
