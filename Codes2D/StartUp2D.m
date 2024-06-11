function S = StartUp2D(S)
% Purpose : Setup script, building operators, grid, metric, and connectivity tables.
% Definition of constants
S.In = 1; S.Out = 2; S.Wall = 3; S.Far = 4; S.Cyl = 5; S.Dirichlet = 6; S.Neuman = 7; S.Slip = 8;
S.Nfp = S.N+1; S.Np = (S.N+1)*(S.N+2)/2; S.Nfaces=3; S.NODETOL = 1e-12;

% Compute nodal set
[S.x,S.y] = Nodes2D(S.N); [S.r,S.s] = xytors(S.x,S.y);

% Build reference element matrices
S.V = Vandermonde2D(S.N,S.r,S.s); S.invV = inv(S.V);
S.MassMatrix = S.invV'*S.invV;
[S.Dr,S.Ds] = Dmatrices2D(S.N, S.r, S.s, S.V);

% build coordinates of all the nodes
va = S.EToV(:,1)'; vb = S.EToV(:,2)'; vc = S.EToV(:,3)';
S.x = 0.5*(-(S.r+S.s)*S.VX(va)+(1+S.r)*S.VX(vb)+(1+S.s)*S.VX(vc));
S.y = 0.5*(-(S.r+S.s)*S.VY(va)+(1+S.r)*S.VY(vb)+(1+S.s)*S.VY(vc));

% find all the nodes that lie on each edge
fmask1   = find( abs(S.s+1) < S.NODETOL)'; 
fmask2   = find( abs(S.r+S.s) < S.NODETOL)';
fmask3   = find( abs(S.r+1) < S.NODETOL)';
S.Fmask  = [fmask1;fmask2;fmask3]';
S.Fx = S.x(S.Fmask(:), :); S.Fy = S.y(S.Fmask(:), :);

% Create surface integral terms
S = Lift2D(S);

% calculate geometric factors
[S.rx,S.sx,S.ry,S.sy,S.J] = GeometricFactors2D(S.x,S.y,S.Dr,S.Ds);

% calculate geometric factors
[S.nx, S.ny, S.sJ] = Normals2D(S);
S.Fscale = S.sJ./(S.J(S.Fmask,:));

% Build connectivity matrix
[S.EToE, S.EToF] = tiConnect2D(S.EToV);

% Build connectivity maps
S=BuildMaps2D(S);

% Compute weak operators (could be done in preprocessing to save time)
[S.Vr, S.Vs] = GradVandermonde2D(S.N, S.r, S.s);
S.Drw = (S.V*S.Vr')/(S.V*S.V'); S.Dsw = (S.V*S.Vs')/(S.V*S.V');


S.rk4a = [            0.0 ...
        -567301805773.0/1357537059087.0 ...
        -2404267990393.0/2016746695238.0 ...
        -3550918686646.0/2091501179385.0  ...
        -1275806237668.0/842570457699.0];
S.rk4b = [ 1432997174477.0/9575080441755.0 ...
         5161836677717.0/13612068292357.0 ...
         1720146321549.0/2090206949498.0  ...
         3134564353537.0/4481467310338.0  ...
         2277821191437.0/14882151754819.0];
S.rk4c = [             0.0  ...
         1432997174477.0/9575080441755.0 ...
         2526269341429.0/6820363962896.0 ...
         2006345519317.0/3224310063776.0 ...
         2802321613138.0/2924317926251.0]; 
S.source=@(x,y,t,U) 0;

return
