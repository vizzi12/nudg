function S=analyticCircleBCtypeNoCurv(S)
% Globals2D;

BCType=zeros(size(S.EToE));

fd = @(p) dcircle(p,0,0,1);
S.BCType = guusCorrectBCTable_v2(S.EToV,S.VX,S.VY,BCType,fd,S.Dirichlet,S.hgrid);


S=BuildBCMaps2D(S);