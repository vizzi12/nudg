function hS=StartE(hL,uL,vL,cL,hR,uR,vR,cR);
g=9.81;
hMin=min(hL,hR);

hS=(1/g)*(0.5*(cL+cR)-0.25*(uR-uL))^2;

if(hS <= hMin)
    return
else
    GEL=sqrt(0.5*g*(hS+hL)/(hS*hL));
    GER=sqrt(0.5*g*(hS+hR)/(hS*hR));
    hS=(GEL*hL+GER*hR-(uR-uL))/(GEL+GER);
end

return;