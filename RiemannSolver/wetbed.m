function [h,u,v]=wetbed(hL,uL,vL,cL,hR,uR,vR,cR,MCELLS,CHALEN,gate,t)
g=9.81;
hS=StartE(hL,uL,vL,cL,hR,uR,vR,cR);
h0=hS;

tol=1e-4;
niter=3000;
j=1;
for i=1:niter
    [FL,FLD]=GeoFun(hS,hL,cL);
    [FR,FRD]=GeoFun(hS,hR,cR);

    hS=hS - (FL+FR+uR-uL)/(FLD+FRD);
%     abs(hS-h0)/(0.5*(hS+h0))
    if (abs(hS-h0)/(0.5*(hS+h0)) < tol)
        j=0;
        break
    end
    if (hS < 0); hS = tol; end
    h0=hS;
end

if (j)
    error('number of iterations exceeded')
end

uS=0.5*(uL+uR) + 0.5*(FR-FL);
cS=sqrt(g*hS);


h=zeros(1,MCELLS);
u=zeros(1,MCELLS);
v=zeros(1,MCELLS);
for i=1:MCELLS
    x=i*CHALEN/MCELLS-gate;
    S=x/t;
    
     [h(i),u(i),v(i)]=SamWet(S,hS,uS,cS,uL,hL,cL,vL,uR,hR,cR,vR);
end
return