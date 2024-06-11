function [p,EToV,VX,VY,K,Nv,hgrid]=guusMesh2D(hgrid)
    if (hgrid == 0)
%         [p,EToV]=load('guusGridh05.mat');
            [p,EToV] = guush05();
            hgrid=0.5;
            [p,EToV]=fixmesh(p,EToV);
    else
        x1=-5;
        x2=15;
        y1=-10;
        y2=10;
        r=1;
        fd = @(p) drectangle(p,x1,x2,y1,y2);
        fd2 = @(p) dcircle(p,0,0,r);
        
        fd = @(p) ddiff(fd(p),fd2(p));
        
        [p,EToV] = distmesh2d(fd,@huniform,hgrid,[x1,y1;x2,y2],[x1,y1;x1,y2;x2,y1;x2,y2]);
        [p,EToV]=fixmesh(p,EToV);
    end
    

    VX = p(:,1)'; % x-coordinates of vertex nodes
    VY = p(:,2)'; % y-coordinates of vertex nodes
    K = size(EToV,1); Nv = length(VX(:));
return