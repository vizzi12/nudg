function [F,FD] = GeoFun(hS,hK,cK)
g=9.81;
if (hS <= hK)
    cS=sqrt(g*hS);
    F=2*(cS-cK);
    FD=g/cS;
else
    GES=sqrt(0.5*g*(hS+hK)/(hS*hK));
    F=(hS-hK)*(GES);
    FD=GES - 0.25*g*(hS-hK)/(GES*hS*hS);
end


return