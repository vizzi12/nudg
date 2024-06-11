function S=analyticBCtype(S)
% Globals2D;

BCType=zeros(size(S.EToE));

x1 = -2; x2 = 2; y1 = -2; y2 = 2;
fd = @(p) drectangle(p,x1,x2,y1,y2);
S.BCType = CorrectBCTable_v2(S.EToV,S.VX,S.VY,BCType,fd,S.Dirichlet);
S=BuildBCMaps2D(S);
