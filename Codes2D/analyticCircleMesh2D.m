function [p,EToV,VX,VY,K,Nv,hgrid]=analyticCircleMesh2D(hgrid)
    if (hgrid == 0)
            [p,EToV] = circle();
            hgrid=0.7;
    else
        fd = @(p) dcircle(p,0,0,1);
        [p,EToV] = distmesh2d(fd,@huniform,hgrid,[-1,-1;1,1],[]);
        [p,EToV]=fixmesh(p,EToV);
    end
    VX = p(:,1)'; % x-coordinates of vertex nodes
    VY = p(:,2)'; % y-coordinates of vertex nodes
    K = size(EToV,1); Nv = length(VX(:));
return