function S=MakeCylinder2D(faces, ra,xo,yo,S)

% Function: MakeCylinder2D(faces, ra, xo, yo)
% Purpose:  Use Gordon-Hall blending with an isoparametric map to modify a list
%           of faces so they conform to a cylinder of radius r centered at (xo,yo)

NCurveFaces = size(faces,1);
vflag = zeros(size(S.VX));
for n=1:NCurveFaces
  % move vertices of faces to be curved onto circle
  k = faces(n,1); f = faces(n,2);
  v1 = S.EToV(k, f); v2 = S.EToV(k, mod(f,S.Nfaces)+1);

  % compute polar angles of start and end face vertices relative to circle center
  theta1 = atan2(S.VY(v1)-yo,S.VX(v1)-xo); 
  theta2 = atan2(S.VY(v2)-yo,S.VX(v2)-xo);

  % move vertices onto circle
  newx1 = xo + ra*cos(theta1); newy1 = yo + ra*sin(theta1);
  newx2 = xo + ra*cos(theta2); newy2 = yo + ra*sin(theta2);

  % update mesh vertex locations
  S.VX(v1) = newx1; S.VX(v2) = newx2; S.VY(v1) = newy1; S.VY(v2) = newy2; 

  % store modified vertex numbers
  vflag(v1) = 1;  vflag(v2) = 1;
end

% map modified vertex flag to each element
vflag = vflag(S.EToV);

% locate elements with at least one modified vertex
ks = find(sum(vflag,2)>0);

% build coordinates of all the corrected nodes
va = S.EToV(ks,1)'; vb = S.EToV(ks,2)'; vc = S.EToV(ks,3)';
S.x(:,ks) = 0.5*(-(S.r+S.s)*S.VX(va)+(1+S.r)*S.VX(vb)+(1+S.s)*S.VX(vc));
S.y(:,ks) = 0.5*(-(S.r+S.s)*S.VY(va)+(1+S.r)*S.VY(vb)+(1+S.s)*S.VY(vc));

for n=1:NCurveFaces  % deform specified faces
  k = faces(n,1); f = faces(n,2);

  % find vertex locations for this face and tangential coordinate
  if(f==1); v1 = S.EToV(k,1); v2 = S.EToV(k,2); vr = S.r; end
  if(f==2); v1 = S.EToV(k,2); v2 = S.EToV(k,3); vr = S.s; end
  if(f==3); v1 = S.EToV(k,1); v2 = S.EToV(k,3); vr = S.s; end
  fr = vr(S.Fmask(:,f));
  x1 = S.VX(v1); y1 = S.VY(v1); x2 = S.VX(v2); y2 = S.VY(v2);

  % move vertices at end points of this face to the cylinder
  theta1 = atan2(y1-yo, x1-xo); theta2 = atan2(y2-yo, x2-xo);

  % check to make sure they are in the same quadrant
  if ((theta2 > 0) && (theta1 < 0) &&  (abs(theta1-theta2) > pi)), theta1 = theta1 + 2*pi; end
  if ((theta1 > 0) && (theta2 < 0) &&  (abs(theta1-theta2) > pi)), theta2 = theta2 + 2*pi; end

  
  % distribute N+1 nodes by arc-length along edge
  theta = 0.5*theta1*(1-fr) + 0.5*theta2*(1+fr);

  % evaluate deformation of coordinates
  fdx = xo + ra*cos(theta)-S.x(S.Fmask(:,f),k); 
  fdy = yo + ra*sin(theta)-S.y(S.Fmask(:,f),k);
  
  % build 1D Vandermonde matrix for face nodes and volume nodes
  Vface = Vandermonde1D(S.N, fr);  Vvol  = Vandermonde1D(S.N, vr);
  % compute unblended volume deformations 
  vdx = Vvol*(Vface\fdx); vdy = Vvol*(Vface\fdy);

  % blend deformation and increment node coordinates
  ids = find(abs(1-vr)>1e-7); % warp and blend
  if(f==1); blend = -(S.r(ids)+S.s(ids))./(1-vr(ids)); end
  if(f==2); blend =      +(S.r(ids)+1)./(1-vr(ids)); end
  if(f==3); blend = -(S.r(ids)+S.s(ids))./(1-vr(ids)); end

  S.x(ids,k) = S.x(ids,k)+blend.*vdx(ids);
  S.y(ids,k) = S.y(ids,k)+blend.*vdy(ids);
end

% repair other coordinate dependent information
S.Fx = S.x(S.Fmask(:), :); S.Fy = S.y(S.Fmask(:), :);
[S.rx,S.sx,S.ry,S.sy,S.J] = GeometricFactors2D(S.x, S.y,S.Dr,S.Ds);
[S.nx, S.ny, S.sJ] = Normals2D(S); S.Fscale = S.sJ./(S.J(S.Fmask,:));
return