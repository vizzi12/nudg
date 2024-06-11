function S=parabolicBCtype(S)
% Globals2D;

BCType=zeros(size(S.EToE));

x1 = -1; x2 = 26; y1 = 0; y2 = 1;
fd = @(p) drectangle(p,x1,x2,y1,y2);
BCType = CorrectBCTable_v2(S.EToV,S.VX,S.VY,BCType,fd,S.Wall);


x1 = 0; x2 = 25; y1 = -1; y2 = 2;
fd = @(p) drectangle(p,x1,x2,y1,y2);
S.BCType = CorrectBCTable_v2(S.EToV,S.VX,S.VY,BCType,fd,S.Dirichlet);

S=BuildBCMaps2D(S);
