function [p,EToV,VX,VY,K,Nv,hgrid]=parabolicMesh2D(hgrid)
    if (hgrid == 0)
        error("Error: No predefined mesh")
    else
        X = 0:hgrid:25;
        Y = 0:hgrid:1;
        [VX,VY]= meshgrid(X,Y); VX = VX(:)'; VY = VY(:)';
        EToV = delaunay(VX,VY);
        K = size(EToV,1); Nv = length(VX(:));
        p=[VX',VX'];
    end

return