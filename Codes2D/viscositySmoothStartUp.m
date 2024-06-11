function S2 = viscositySmoothStartUp(S)

S2.N=2;
S2.p=S.p; S2.EToV=S.EToV; S2.VX=S.VX; S2.VY=S.VY; S2.K=S.K; S2.Nv=S.Nv; S2.hgrid=S.hgrid;
S2=StartUp2D(S2);

S2.counts=groupcounts(S.EToV(:));
S2.id=repelem(1:S.Nv,S2.counts)';
S2.ind=[];
for i=1:S.Nv
    S2.ind=[S2.ind;find( (S2.x-S2.VX(i)).^2+(S2.y-S2.VY(i)).^2 < S2.NODETOL  )];
end

return