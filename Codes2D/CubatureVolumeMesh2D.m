function cub = CubatureVolumeMesh2D(CubatureOrder,S)

% function cub = CubatureVolumeMesh2D(CubatureOrder)
% purpose: build cubature nodes, weights and geometric factors for all elements

% Globals2D;

% set up cubature nodes
[cub.R,cub.S,cub.W, cub.Ncub] = Cubature2D(CubatureOrder); 



% evaluate generalized Vandermonde of Lagrange interpolant functions at cubature nodes
cub.V  = InterpMatrix2D(cub.R, cub.S,S); 

% evaluate local derivatives of Lagrange interpolants at cubature nodes
[cub.Dr,cub.Ds] = Dmatrices2D(S.N,cub.R,cub.S,S.V);

% evaluate the geometric factors at the cubature nodes
[cub.rx,cub.sx,cub.ry,cub.sy,cub.J] = GeometricFactors2D(S.x,S.y, cub.Dr,cub.Ds);

% custom mass matrix per element
cub.mmCHOL = zeros(S.Np, S.Np, S.K); cub.mm = zeros(S.Np, S.Np, S.K);
for k=1:S.K
  cub.mm(:,:,k)     = cub.V'*diag(cub.J(:,k).*cub.W(:))*cub.V;
  cub.mmCHOL(:,:,k) = chol(cub.mm(:,:,k));
end

% incorporate weights and Jacobian
cub.w = cub.W; cub.W = cub.W*ones(1,S.K); cub.W = cub.W.*cub.J; 

% compute coordinates of cubature nodes
cub.x = cub.V*S.x; cub.y = cub.V*S.y;
return
