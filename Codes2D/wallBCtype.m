function S=wallBCtype(S)
% Globals2D;

BCType=zeros(size(S.EToE));

x1 = 0; x2 = 200; y1 = -100; y2 = 300;
fd = @(p) drectangle(p,x1,x2,y1,y2);
BCType = CorrectBCTable_v2(S.EToV,S.VX,S.VY,BCType,fd,S.Dirichlet);

x1 = -100; x2 = 300; y1 = 0; y2 = 200;
fd = @(p) drectangle(p,x1,x2,y1,y2);
fd2 = @(p) drectangle(p,95,105,0,87.5);
fd3 = @(p) drectangle(p,95,105,162.5,200);

fd = @(p) ddiff(fd(p),fd2(p));
fd = @(p) ddiff(fd(p),fd3(p));

S.BCType = CorrectBCTable_v2(S.EToV,S.VX,S.VY,BCType,fd,S.Wall);
S=BuildBCMaps2D(S);