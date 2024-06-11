function S=circularBCtype(S)

S.BCType = zeros(size(S.EToE));

x1 = -20; x2 = 20; y1 = -20; y2 = 20;
fd = @(p) drectangle(p,x1,x2,y1,y2);
S.BCType = CorrectBCTable_v2(S.EToV,S.VX,S.VY,S.BCType,fd,S.Dirichlet);
S=BuildBCMaps2D(S);
return;