function [LIFT] = Lift2D(G)
Emat = zeros(G.Np, G.Nfaces*G.Nfp);

% face 1
faceR = G.r(G.Fmask(:,1));
V1D = Vandermonde1D(G.N, faceR); 
massEdge1 = inv(V1D*V1D');
Emat(G.Fmask(:,1),1:G.Nfp) = massEdge1;

% face 2
faceR = G.r(G.Fmask(:,2));
V1D = Vandermonde1D(G.N, faceR);
massEdge2 = inv(V1D*V1D');
Emat(G.Fmask(:,2),G.Nfp+1:2*G.Nfp) = massEdge2;

% face 3
faceS = G.s(G.Fmask(:,3));
V1D = Vandermonde1D(G.N, faceS); 
massEdge3 = inv(V1D*V1D');
Emat(G.Fmask(:,3),2*G.Nfp+1:3*G.Nfp) = massEdge3;

% inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
LIFT = G.V*(G.V'*Emat);
return
