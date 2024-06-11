function [p,EToV,VX,VY,K,Nv,hgrid]=circularMesh2D(hgrid)
    if (hgrid == 0)
        load('DamBreakGridh02'); 
        hgrid=0.2;
    else
        fd = @(p) drectangle(p,-20,20,-20,20);
        [p,EToV] = distmesh2d(fd,@huniform,hgrid,[-20,-20;20,20],[-20,-20;-20,20;20,-20;20,20]);
        [p,EToV]=fixmesh(p,EToV);
    end
    VX = p(:,1)'; % x-coordinates of vertex nodes
    VY = p(:,2)'; % y-coordinates of vertex nodes
    K = size(EToV,1); Nv = length(VX(:));
return