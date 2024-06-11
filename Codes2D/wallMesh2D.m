function [p,EToV,VX,VY,K,Nv,hgrid]=wallMesh2D(hgrid)
    if (hgrid == 0)
        load("DamBreakEskilssonGridh5.mat"); 
        hgrid=5;
%     [p,EToV] = wall();
%     hgrid=5;
    else
        fd = @(p) drectangle(p,0,200,0,200);
        fd2 = @(p) drectangle(p,95,105,0,87.5);
        fd3 = @(p) drectangle(p,95,105,162.5,200);
        
        fd = @(p) ddiff(fd(p),fd2(p));
        fd = @(p) ddiff(fd(p),fd3(p));
        
        yv = [0:hgrid:87.5,87.5:hgrid:162.5,162.5:hgrid:200]';
        xv = 0*yv+105;
        wPoints = [xv,yv];
    %     fh = @(p) min( (abs(p(:,1)-100)+40)/100 , 1 );
        [p,EToV] = distmesh2d(fd,@huniform,hgrid,[0,0;200,200],[[0,0;0,200;200,0;200,200;95,0;105,0;95,87.5;105,87.5;95,162.5;105,162.5;95,200;105,200];wPoints]);
        [p,EToV]=fixmesh(p,EToV);
    end
    

    VX = p(:,1)'; % x-coordinates of vertex nodes
    VY = p(:,2)'; % y-coordinates of vertex nodes
    K = size(EToV,1); Nv = length(VX(:));
return