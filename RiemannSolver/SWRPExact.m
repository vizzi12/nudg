function [h,u,v]=SWRPExact(hL,uL,vL,hR,uR,vR,MCELLS,CHALEN,gate,t)
g=9.81;

cL=sqrt(g*hL);
cR=sqrt(g*hR);

DCrit=(uR-uL) - 2*(cL+cR);

if( (hL <= 0) || (hR <= 0) || (DCrit >= 0) )
    error('Dry bed')
else
    [h,u,v]=wetbed(hL,uL,vL,cL,hR,uR,vR,cR,MCELLS,CHALEN,gate,t);
end


return