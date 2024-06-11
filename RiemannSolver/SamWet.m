function [h,u,v]=SamWet(S,hS,uS,cS,uL,hL,cL,vL,uR,hR,cR,vR)
g=9.81;
if(S <= uS)
%Left wave
    if (hS >= hL)
        QL=sqrt((hS+hL)*hS/(2*hL*hL));
        SL=uL-cL*QL;
        if(S <= SL)
            h = hL;
            u = uL;
        else
            h = hS;
            u = uS; 
        end
    else
        SHL=uL-cL;
        if(S <= SHL)
            h=hL;
            u=uL;
        else
            STL= uS-cS;
            if(S <= STL)
                u=(uL+2*cL+2*S)/3;
                C=(uL+2*cL-S)/3;
                h=C*C/g;
            else
                h=hS;
                u=uS;
            end
        end
    end
    v=vL;
else
%Right wave
    if (hS >= hR)
        QR=sqrt((hS+hR)*hS/(2*hR*hR));
        SR=uR+cR*QR;
        if(S >= SR)
            h = hR;
            u = uR;
        else
            h = hS;
            u = uS; 
        end
    else
        SHR=uR+cR;
        if(S >= SHR)
            h=hR;
            u=uR;
        else
            STR= uS+cS;
            if(S >= STR)
                u=(uR-2*cR+2*S)/3;
                C=(-uR+2*cR+S)/3;
                h=C*C/g;
            else
                h=hS;
                u=uS;
            end
        end
    end
    v=vR;
end





return