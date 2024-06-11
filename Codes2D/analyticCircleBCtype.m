function S=analyticCircleBCtype(S)
% Globals2D;

BCType=zeros(size(S.EToE));

fd = @(p) dcircle(p,0,0,1);
S.BCType = guusCorrectBCTable_v2(S.EToV,S.VX,S.VY,BCType,fd,S.Cyl,S.hgrid);

[k,f] = find(S.BCType==S.Cyl);
S.curved = [];
if(~isempty(k))
    cylfaces = [k,f];
    S.curved = sort(unique(k)); 
    S=MakeCylinder2D(cylfaces, 1, 0, 0,S);
    %turn cylinders into walls
    ids = find(S.BCType==S.Cyl);  S.BCType(ids) = S.Dirichlet;
end
S.straight = setdiff(1:S.K, S.curved);

S=BuildBCMaps2D(S);