function [p,EToV,VX,VY,K,Nv,hgrid]=analyticMesh2D(hgrid)
    if (hgrid == 0)
        error("Error: No predefined mesh")
    else
        fd = @(p) drectangle(p,-2,2,-2,2);
        [p,EToV] = distmesh2d(fd,@huniform,hgrid,[-2,-2;2,2],[-2,-2;-2,2;2,-2;2,2]);
        [p,EToV]=fixmesh(p,EToV);
    end
    VX = p(:,1)'; % x-coordinates of vertex nodes
    VY = p(:,2)'; % y-coordinates of vertex nodes
    K = size(EToV,1); Nv = length(VX(:));
return