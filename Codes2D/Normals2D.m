function [nx, ny, sJ] = Normals2D(S)

% function [nx, ny, sJ] = Normals2D()
% Purpose : Compute outward pointing normals at elements faces and surface Jacobians

% Globals2D;
xr = S.Dr*S.x; yr = S.Dr*S.y; xs = S.Ds*S.x; ys = S.Ds*S.y; %J = S.xr.*ys-xs.*yr;

% interpolate geometric factors to face nodes
fxr = xr(S.Fmask, :); fxs = xs(S.Fmask, :); fyr = yr(S.Fmask, :); fys = ys(S.Fmask, :);

% build normals
nx = zeros(3*S.Nfp, S.K); ny = zeros(3*S.Nfp, S.K);
fid1 = (1:S.Nfp)'; fid2 = (S.Nfp+1:2*S.Nfp)'; fid3 = (2*S.Nfp+1:3*S.Nfp)';

% face 1
nx(fid1, :) =  fyr(fid1, :); ny(fid1, :) = -fxr(fid1, :);

% face 2
nx(fid2, :) =  fys(fid2, :)-fyr(fid2, :); ny(fid2, :) = -fxs(fid2, :)+fxr(fid2, :);

% face 3
nx(fid3, :) = -fys(fid3, :); ny(fid3, :) =  fxs(fid3, :);

% normalise
sJ = sqrt(nx.*nx+ny.*ny); nx = nx./sJ; ny = ny./sJ;
return;
