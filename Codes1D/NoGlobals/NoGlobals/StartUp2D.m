function [G] = StartUp2D(G)
G.Nfp = G.N+1; G.Np = (G.N+1)*(G.N+2)/2; G.Nfaces=3; G.NODETOL = 1e-12;

% Compute nodal set
[G.x,G.y] = Nodes2D(G.N); [G.r,G.s] = xytors(G.x,G.y);

% Build reference element matrices
G.V = Vandermonde2D(G.N,G.r,G.s); invV = inv(G.V);
G.MassMatrix = invV'*invV;
[G.Dr,G.Ds] = Dmatrices2D(G.N, G.r, G.s, G.V);

% build coordinates of all the nodes
va = G.EToV(:,1)'; vb = G.EToV(:,2)'; vc = G.EToV(:,3)';
G.x = 0.5*(-(G.r+G.s)*G.VX(va)+(1+G.r)*G.VX(vb)+(1+G.s)*G.VX(vc));
G.y = 0.5*(-(G.r+G.s)*G.VY(va)+(1+G.r)*G.VY(vb)+(1+G.s)*G.VY(vc));

% find all the nodes that lie on each edge
fmask1   = find( abs(G.s+1) < G.NODETOL)'; 
fmask2   = find( abs(G.r+G.s) < G.NODETOL)';
fmask3   = find( abs(G.r+1) < G.NODETOL)';
G.Fmask  = [fmask1;fmask2;fmask3]';
G.Fx = G.x(G.Fmask(:), :); G.Fy = G.y(G.Fmask(:), :);

% Create surface integral terms
G.LIFT = Lift2D(G);

% calculate geometric factors
[G.rx,G.sx,G.ry,G.sy,G.J] = GeometricFactors2D(G.x,G.y,G.Dr,G.Ds);

% calculate geometric factors
[G.nx, G.ny, G.sJ] = Normals2D(G);
G.Fscale = G.sJ./(G.J(G.Fmask,:));

% Build connectivity matrix
[G.EToE, G.EToF] = tiConnect2D(G.EToV);

% Build connectivity maps
[G.mapM, G.mapP, G.vmapM, G.vmapP, G.vmapB, G.mapB] = BuildMaps2D(G);

% Compute weak operators (could be done in preprocessing to save time)
[Vr, Vs] = GradVandermonde2D(G.N, G.r, G.s);
G.Drw = (G.V*Vr')/(G.V*G.V'); G.Dsw = (G.V*Vs')/(G.V*G.V');

G.Dx = G.Dr*G.rx + G.Ds*G.sx;
G.Dy = G.Dr*G.ry + G.Ds*G.sy;

G.rk4a = [            0.0 ...
        -567301805773.0/1357537059087.0 ...
        -2404267990393.0/2016746695238.0 ...
        -3550918686646.0/2091501179385.0  ...
        -1275806237668.0/842570457699.0];
G.rk4b = [ 1432997174477.0/9575080441755.0 ...
         5161836677717.0/13612068292357.0 ...
         1720146321549.0/2090206949498.0  ...
         3134564353537.0/4481467310338.0  ...
         2277821191437.0/14882151754819.0];
G.rk4c = [             0.0  ...
         1432997174477.0/9575080441755.0 ...
         2526269341429.0/6820363962896.0 ...
         2006345519317.0/3224310063776.0 ...
         2802321613138.0/2924317926251.0]; 

end

