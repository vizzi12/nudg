function [mapM, mapP, vmapM, vmapP, vmapB, mapB] = BuildMaps2D(G)

% function [mapM, mapP, vmapM, vmapP, vmapB, mapB] = BuildMaps2D
% Purpose: Connectivity and boundary tables in the K # of Np elements

% number volume nodes consecutively
nodeids = reshape(1:G.K*G.Np, G.Np, G.K);
vmapM   = zeros(G.Nfp, G.Nfaces, G.K); vmapP   = zeros(G.Nfp, G.Nfaces, G.K); 
mapM    = (1:G.K*G.Nfp*G.Nfaces)';     mapP    = reshape(mapM, G.Nfp, G.Nfaces, G.K);
 
% find index of face nodes with respect to volume node ordering
for k1=1:G.K
  for f1=1:G.Nfaces
    vmapM(:,f1,k1) = nodeids(G.Fmask(:,f1), k1);
  end
end

one = ones(1, G.Nfp);
for k1=1:G.K
  for f1=1:G.Nfaces
    % find neighbor
    k2 = G.EToE(k1,f1); f2 = G.EToF(k1,f1);
    
    % reference length of edge
    v1 = G.EToV(k1,f1); v2 = G.EToV(k1, 1+mod(f1,G.Nfaces));
    refd = sqrt( (G.VX(v1)-G.VX(v2))^2 + (G.VY(v1)-G.VY(v2))^2 );

    % find find volume node numbers of left and right nodes 
    vidM = vmapM(:,f1,k1); vidP = vmapM(:,f2,k2);    
    x1 = G.x(vidM); y1 = G.y(vidM); x2 = G.x(vidP); y2 = G.y(vidP);
    x1 = x1*one;  y1 = y1*one;  x2 = x2*one;  y2 = y2*one;

    % Compute distance matrix
    D = (x1 -x2').^2 + (y1-y2').^2;
    [idM, idP] = find(sqrt(abs(D))<G.NODETOL*refd);
    vmapP(idM,f1,k1) = vidP(idP); mapP(idM,f1,k1) = idP + (f2-1)*G.Nfp+(k2-1)*G.Nfaces*G.Nfp;
  end
end

% reshape vmapM and vmapP to be vectors and create boundary node list
vmapP = vmapP(:); vmapM = vmapM(:); mapP = mapP(:); 
mapB = find(vmapP==vmapM); vmapB = vmapM(mapB);
return
