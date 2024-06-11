function  gauss = GaussFaceMesh2D(NGauss,S)

% function: gauss = GaussFaceMesh2D(NGauss)
% Purpose:  Compute Gauss nodes for face term integration, and interpolation matrices 

% Globals2D;
gauss.NGauss = NGauss;

[gauss.z, gauss.w] = JacobiGQ(0, 0, NGauss-1);
face1r =  gauss.z;         face2r = -gauss.z; face3r = -ones(NGauss, 1); 
face1s = -ones(NGauss, 1); face2s =  gauss.z; face3s = -gauss.z;

gauss.finterp = zeros(NGauss,S.Np,S.Nfaces);
V1 = Vandermonde2D(S.N, face1r, face1s);   gauss.finterp(:,:,1) = V1*S.invV;
V2 = Vandermonde2D(S.N, face2r, face2s);   gauss.finterp(:,:,2) = V2*S.invV;
V3 = Vandermonde2D(S.N, face3r, face3s);   gauss.finterp(:,:,3) = V3*S.invV;

gauss.interp = [gauss.finterp(:,:,1); gauss.finterp(:,:,2); gauss.finterp(:,:,3)];
gauss.mapM = reshape(1:NGauss*S.Nfaces*S.K, NGauss*S.Nfaces, S.K);		    
gauss.mapP = reshape(1:NGauss*S.Nfaces*S.K, NGauss*S.Nfaces, S.K);		    
gauss.mapI = []; gauss.mapO = []; gauss.mapW = []; gauss.mapB = [];
gauss.mapS = []; gauss.mapD = []; gauss.mapN = []; gauss.mapC = [];

zer = zeros(NGauss*S.Nfaces,S.K);
gauss.nx = zer; gauss.ny = zer; gauss.rx = zer; gauss.ry = zer; gauss.sx = zer; gauss.sy = zer; gauss.J = zer; gauss.sJ = zer;
for f1=1:S.Nfaces
  VM = gauss.finterp(:,:,f1);
  dVMdr = VM*S.Dr;  dVMds = VM*S.Ds;
  ids1 = (f1-1)*NGauss+1:f1*NGauss;
  for k1=1:S.K
    % calculate geometric factors at Gauss points
    [grx,gsx,gry,gsy,gJ] = ...
	GeometricFactors2D(S.x(:,k1),S.y(:,k1),dVMdr,dVMds);
    % compute normals at Gauss points
    if(f1==1)  gnx = -gsx;       gny = -gsy;      end;
    if(f1==2)  gnx =  grx+gsx;   gny =  gry+gsy;  end;
    if(f1==3)  gnx = -grx;       gny = -gry;      end;
    
    gsJ = sqrt(gnx.*gnx+gny.*gny);
    gnx = gnx./gsJ;  gny = gny./gsJ;  gsJ = gsJ.*gJ;
    
    gauss.nx(ids1,k1) = gnx; gauss.ny(ids1,k1) = gny; gauss.sJ(ids1,k1) = gsJ;
    gauss.rx(ids1,k1) = grx; gauss.ry(ids1,k1) = gry; gauss.J (ids1,k1) = gJ;
    gauss.sx(ids1,k1) = gsx; gauss.sy(ids1,k1) = gsy;

    k2 = S.EToE(k1,f1); f2 = S.EToF(k1,f1); ids2 = f2*NGauss:-1:(f2-1)*NGauss+1;
    if(k1~=k2)
      gauss.mapP(ids1, k1) = gauss.mapM(ids2,k2);		   
    else
      gauss.mapP(ids1, k1) = gauss.mapM(ids1,k1);	
      gauss.mapB = [gauss.mapB;gauss.mapM(ids1,k1)]; 
      switch(S.BCType(k1,f1))
	case {S.Wall}
	  gauss.mapW = [gauss.mapW;gauss.mapM(ids1,k1)]; 
	case {S.In}
	  gauss.mapI = [gauss.mapI;gauss.mapM(ids1,k1)]; 
	case {S.Out}
	  gauss.mapO = [gauss.mapO;gauss.mapM(ids1,k1)]; 
	case {S.Slip}
	  gauss.mapS = [gauss.mapS;gauss.mapM(ids1,k1)]; 
	case {S.Dirichlet}
	  gauss.mapD = [gauss.mapD;gauss.mapM(ids1,k1)]; 
	case {S.Neuman}
	  gauss.mapN = [gauss.mapN;gauss.mapM(ids1,k1)]; 
	case {S.Cyl}
	  gauss.mapC = [gauss.mapC;gauss.mapM(ids1,k1)]; 
      end
    end
  end
end
gauss.x = gauss.interp*S.x;
gauss.y = gauss.interp*S.y;
gauss.W = [gauss.w;gauss.w;gauss.w]*ones(1,S.K);
gauss.W = gauss.W.*gauss.sJ;

gauss.mapM3D=cat(3,gauss.mapM,gauss.mapM+length(gauss.mapM(:)),gauss.mapM+2*length(gauss.mapM(:)));
gauss.mapP3D=cat(3,gauss.mapP,gauss.mapP+length(gauss.mapP(:)),gauss.mapP+2*length(gauss.mapP(:)));

return;
