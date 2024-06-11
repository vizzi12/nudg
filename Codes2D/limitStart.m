function limit = limitStart(S)
% Globals2D
AVE = sum(S.MassMatrix)/2;

% Compute displacements from center of nodes for Taylor expansion of limited fields

% Find neighbors in patch
E1 = S.EToE(:,1)'; E2 = S.EToE(:,2)'; E3 = S.EToE(:,3)';

% Extract coordinates of vertices and centers of elements
v1 = S.EToV(:,1); xv1 = S.VX(v1); yv1 = S.VY(v1);
v2 = S.EToV(:,2); xv2 = S.VX(v2); yv2 = S.VY(v2);
v3 = S.EToV(:,3); xv3 = S.VX(v3); yv3 = S.VY(v3);  

% compute face unit normals and lengths
fnx = [yv2-yv1;yv3-yv2;yv1-yv3]; fny = -[xv2-xv1;xv3-xv2;xv1-xv3];
fL = sqrt(fnx.^2 + fny.^2); fnx = fnx./fL; fny = fny./fL;

% compute element centers
xc0 = AVE*S.x; xc1 = xc0(E1); xc2 = xc0(E2); xc3 = xc0(E3);
yc0 = AVE*S.y; yc1 = yc0(E1); yc2 = yc0(E2); yc3 = yc0(E3);

% Compute weights for face gradients 
A0 = AVE*S.J*2/3;

% Find boundary faces for each face 
id1 = find(S.BCType(:,1)); id2 = find(S.BCType(:,2)); id3 = find(S.BCType(:,3));

% Compute location of centers of reflected ghost elements at boundary faces
H1 = 2*(A0(id1)./fL(1,id1)); 
xc1(id1) = xc1(id1) + 2*fnx(1,id1).*H1; yc1(id1) = yc1(id1) + 2*fny(1,id1).*H1; 

H2 = 2*(A0(id2)./fL(2,id2)); 
xc2(id2) = xc2(id2) + 2*fnx(2,id2).*H2; yc2(id2) = yc2(id2) + 2*fny(2,id2).*H2; 

H3 = 2*(A0(id3)./fL(3,id3)); 
xc3(id3) = xc3(id3) + 2*fnx(3,id3).*H3; yc3(id3) = yc3(id3) + 2*fny(3,id3).*H3; 

if mod(S.Nfp,2) 
    ids=[ceil(S.Nfp/2),S.Nfp+ceil(S.Nfp/2),2*S.Nfp+ceil(S.Nfp/2)];
    vmapM1 = reshape(S.vmapM, S.Nfp*S.Nfaces, S.K); vmapM1 = vmapM1(ids, :);
    mx=S.x(vmapM1);
    my=S.y(vmapM1);
else
    ids=[S.Nfp/2,S.Nfp/2+1,S.Nfp+S.Nfp/2,S.Nfp+S.Nfp/2+1,2*S.Nfp+S.Nfp/2,2*S.Nfp+S.Nfp/2+1];
    vmapM1 = reshape(S.vmapM, S.Nfp*S.Nfaces, S.K); vmapM1 = vmapM1(ids, :);
    mx=( S.x(vmapM1(1:2:5,:))+S.x(vmapM1(2:2:6,:)) ) ./ 2 ;
    my=( S.y(vmapM1(1:2:5,:))+S.y(vmapM1(2:2:6,:)) ) ./ 2 ;
end

mx=reshape(mx,3,1,[]);
my=reshape(my,3,1,[]);

xc0=reshape(xc0,1,1,[]); xc1=reshape(xc1,1,1,[]); xc2=reshape(xc2,1,1,[]); xc3=reshape(xc3,1,1,[]);
yc0=reshape(yc0,1,1,[]); yc1=reshape(yc1,1,1,[]); yc2=reshape(yc2,1,1,[]); yc3=reshape(yc3,1,1,[]);

alpha1=pagemldivide([xc1-xc0,xc2-xc0 ; yc1-yc0,yc2-yc0],[mx(1,:,:)-xc0;my(1,:,:)-yc0]);
alpha2=pagemldivide([xc2-xc0,xc3-xc0 ; yc2-yc0,yc3-yc0],[mx(2,:,:)-xc0;my(2,:,:)-yc0]);
alpha3=pagemldivide([xc3-xc0,xc1-xc0 ; yc3-yc0,yc1-yc0],[mx(3,:,:)-xc0;my(3,:,:)-yc0]);


mx=reshape(mx,1,3,[]);
my=reshape(my,1,3,[]);

limit.s=[mx-xc0;my-yc0];

limit.s=limit.s./vecnorm(limit.s);
limit.s=reshape(permute(limit.s,[1,3,2]),[1,2,S.K,3]);

limit.mx=mx;
limit.my=my;

limit.xc0=xc0; limit.xc1=xc1; limit.xc2=xc2; limit.xc3=xc3;
limit.yc0=yc0; limit.yc1=yc1; limit.yc2=yc2; limit.yc3=yc3;

limit.alpha=permute(cat(2,alpha1,alpha2,alpha3),[2,1,3]);

limit.AVE=AVE;

limit.fnx=fnx;
limit.fny=fny;

limit.idW = find(S.BCType'==S.Wall); limit.idD = find(S.BCType'==S.Dirichlet); 

return