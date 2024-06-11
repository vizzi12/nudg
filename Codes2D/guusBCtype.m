function S=guusBCtype(S)
% Globals2D;

% x1=-10;
% x2=10;
% y1=-5;
% y2=5;
% r=2;

BCType=zeros(size(S.EToE));

x1 = -5; x2 = 16; y1 = -11; y2 = 11;
fd = @(p) drectangle(p,x1,x2,y1,y2);
BCType = CorrectBCTable_v2(S.EToV,S.VX,S.VY,BCType,fd,S.Dirichlet);

% x1 = -6; x2 = 15; y1 = -11; y2 = 11;
% fd = @(p) drectangle(p,x1,x2,y1,y2);
% BCType = CorrectBCTable_v2(S.EToV,S.VX,S.VY,BCType,fd,S.Out);

x1 = -6; x2 = 16; y1 = -10; y2 = 10;
fd = @(p) drectangle(p,x1,x2,y1,y2);

BCType = CorrectBCTable_v2(S.EToV,S.VX,S.VY,BCType,fd,S.Wall);

fd = @(p) dcircle(p,0,0,1);
S.BCType = guusCorrectBCTable_v2(S.EToV,S.VX,S.VY,BCType,fd,S.Cyl,S.hgrid);

        [k,f] = find(S.BCType==S.Cyl);
        S.curved = [];
        if(~isempty(k))
            cylfaces = [k,f];
            S.curved = sort(unique(k)); 
            S=MakeCylinder2D(cylfaces, 1, 0, 0,S);
            %turn cylinders into walls
            ids = find(S.BCType==S.Cyl);  S.BCType(ids) = S.Wall;
        end
        S.straight = setdiff(1:S.K, S.curved);

% S.straight=1:S.K;
% S.curved=[];


S=BuildBCMaps2D(S);
