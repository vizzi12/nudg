function S = BuildMaps2D(S)

% function [mapM, mapP, vmapM, vmapP, vmapB, mapB] = BuildMaps2D
% Purpose: Connectivity and boundary tables in the K # of Np elements

% Globals2D;

% number volume nodes consecutively
nodeids = reshape(1:S.K*S.Np, S.Np, S.K);
vmapM   = zeros(S.Nfp, S.Nfaces, S.K); vmapP   = zeros(S.Nfp, S.Nfaces, S.K); 
mapM    = (1:S.K*S.Nfp*S.Nfaces)';     mapM    = reshape(mapM, S.Nfp, S.Nfaces, S.K);
mapP=zeros(size(mapM));
% find index of face nodes with respect to volume node ordering
for k1=1:S.K
  for f1=1:S.Nfaces
    vmapM(:,f1,k1) = nodeids(S.Fmask(:,f1), k1);
  end
end

one = ones(1, S.Nfp);
for k1=1:S.K
  for f1=1:S.Nfaces
    % find neighbor
    k2 = S.EToE(k1,f1); f2 = S.EToF(k1,f1);
    
    % reference length of edge
    v1 = S.EToV(k1,f1); v2 = S.EToV(k1, 1+mod(f1,S.Nfaces));
    refd = sqrt( (S.VX(v1)-S.VX(v2))^2 + (S.VY(v1)-S.VY(v2))^2 );

    % find find volume node numbers of left and right nodes 
    vidM = vmapM(:,f1,k1); vidP = vmapM(:,f2,k2);    
    x1 = S.x(vidM); y1 = S.y(vidM); x2 = S.x(vidP); y2 = S.y(vidP);
    x1 = x1*one;  y1 = y1*one;  x2 = x2*one;  y2 = y2*one;

    % Compute distance matrix
    D = (x1 -x2').^2 + (y1-y2').^2;
    [idM, idP] = find(sqrt(abs(D))<S.NODETOL*refd);
    vmapP(idM,f1,k1) = vidP(idP); mapP(idM,f1,k1) = idP + (f2-1)*S.Nfp+(k2-1)*S.Nfaces*S.Nfp;
  end
end

% reshape vmapM and vmapP to be vectors and create boundary node list
S.vmapP = vmapP(:); S.vmapM = vmapM(:); S.mapP = mapP(:); 
S.vmapP = vmapP; S.vmapM = vmapM; S.mapP = mapP; S.mapM = mapM;  
S.mapB = find(S.vmapP==S.vmapM); S.vmapB = vmapM(S.mapB);

S.vmapM = reshape(S.vmapM, S.Nfp*S.Nfaces, S.K); S.vmapP = reshape(S.vmapP, S.Nfp*S.Nfaces, S.K);
return
