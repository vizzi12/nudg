function [nx, ny, sJ] = Normals2D(G)

% function [nx, ny, sJ] = Normals2D()
% Purpose : Compute outward pointing normals at elements faces and surface Jacobians

xr = G.Dr*G.x; yr = G.Dr*G.y; xs = G.Ds*G.x; ys = G.Ds*G.y; J = xr.*ys-xs.*yr;

% interpolate geometric factors to face nodes
fxr = xr(G.Fmask, :); fxs = xs(G.Fmask, :); fyr = yr(G.Fmask, :); fys = ys(G.Fmask, :);

% build normals
nx = zeros(3*G.Nfp, G.K); ny = zeros(3*G.Nfp, G.K);
fid1 = (1:G.Nfp)'; fid2 = (G.Nfp+1:2*G.Nfp)'; fid3 = (2*G.Nfp+1:3*G.Nfp)';

% face 1
nx(fid1, :) =  fyr(fid1, :); ny(fid1, :) = -fxr(fid1, :);

% face 2
nx(fid2, :) =  fys(fid2, :)-fyr(fid2, :); ny(fid2, :) = -fxs(fid2, :)+fxr(fid2, :);

% face 3
nx(fid3, :) = -fys(fid3, :); ny(fid3, :) =  fxs(fid3, :);

% normalise
sJ = sqrt(nx.*nx+ny.*ny); nx = nx./sJ; ny = ny./sJ;
return;
