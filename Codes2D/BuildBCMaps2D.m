function [S]=BuildBCMaps2D(S)

% function BuildMaps2DBC
% Purpose: Build specialized nodal maps for various types of
%          boundary conditions, specified in BCType. 

% Globals2D;

% create label of face nodes with boundary types from BCType
bct    = S.BCType';
bnodes = ones(S.Nfp, 1)*bct(:)';
bnodes = bnodes(:);

% find location of boundary nodes in face and volume node lists
S.mapI = find(bnodes==S.In);           S.vmapI = S.vmapM(S.mapI);
S.mapO = find(bnodes==S.Out);          S.vmapO = S.vmapM(S.mapO);
S.mapW = find(bnodes==S.Wall);         S.vmapW = S.vmapM(S.mapW);
S.mapF = find(bnodes==S.Far);          S.vmapF = S.vmapM(S.mapF);
S.mapC = find(bnodes==S.Cyl);          S.vmapC = S.vmapM(S.mapC);
S.mapD = find(bnodes==S.Dirichlet);    S.vmapD = S.vmapM(S.mapD);
S.mapN = find(bnodes==S.Neuman);       S.vmapN = S.vmapM(S.mapN);
S.mapS = find(bnodes==S.Slip);         S.vmapS = S.vmapM(S.mapS);
return;
