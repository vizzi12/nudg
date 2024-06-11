function flux = SWEHLLC2D(nx, ny, UM, UP,form)
g=9.81;

hM = UM(:,:,1); qM = UM(:,:,2); pM = UM(:,:,3);
hP = UP(:,:,1); qP = UP(:,:,2); pP = UP(:,:,3);

qMr = nx.*qM + ny.*pM;
pMr = -ny.*qM + nx.*pM;

qPr = nx.*qP + ny.*pP;
pPr = -ny.*qP + nx.*pP;

aM=sqrt(g*hM);
aP=sqrt(g*hP);

uMr=qMr./hM;
vMr=pMr./hM;

uPr=qPr./hP;
vPr=pPr./hP;


hS=1/g*(0.5*(aM+aP) + 0.25 *(uMr-uPr)).^2;

qL=sqrt( 0.5*(hS+hM).*hS./hM.^2 );
qL(hS <= hM) = 1;

qR=sqrt( 0.5*(hS+hP).*hS./hP.^2 );
qR(hS <= hP) = 1;

sL=uMr-aM.*qL;
sR=uPr+aP.*qR;

sS=( sL.*hP.*(uPr-sR) - sR.*hM.*(uMr-sL) ) ./ ( hP.*(uPr-sR) - hM.*(uMr-sL));

FLh= qMr;
FRh= qPr;

QLh= hM.*( (sL-uMr) ./ (sL - sS)  );
QRh= hP.*( (sR-uPr) ./ (sR - sS)  );

idx1= sS >= 0 & sL < 0 & sR > 0;
idx2= sS < 0 & sL < 0 & sR > 0;
idx3= sL >= 0;
idx4= sR <= 0;

fh=idx1.*( FLh+sL.*(QLh-hM) ) + idx2.*(FRh+sR.*(QRh-hP) )+idx3.*FLh+idx4.*FRh;

FLq= qMr.*uMr+0.5*g*hM.^2;
FRq= qPr.*uPr+0.5*g*hP.^2;

QLq= QLh.*sS;
QRq= QRh.*sS;

fq=idx1.*( FLq+sL.*(QLq-qMr) ) + idx2.*(FRq+sR.*(QRq-qPr) )+idx3.*FLq+idx4.*FRq;

FLp= uMr.*vMr.*hM;
FRp= uPr.*vPr.*hP;

QLp= QLh.*vMr;
QRp= QRh.*vPr;

fp=idx1.*( FLp+sL.*(QLp-pMr) ) + idx2.*(FRp+sR.*(QRp-pPr) )+idx3.*FLp+idx4.*FRp;

flux=zeros(size(UM));
switch form
    case {'weak'}
        flux(:,:,1) = fh;
        flux(:,:,2) = nx.*fq - ny.*fp;
        flux(:,:,3) = ny.*fq + nx.*fp;
    case {'strong'}
        fhd=qMr-fh;
        fqd=qMr.^2./hM+0.5*g*hM.^2 - fq;
        fpd=(qMr.*pMr./hM) - fp;

        flux(:,:,1) = fhd;
        flux(:,:,2) = nx.*fqd - ny.*fpd;
        flux(:,:,3)  =ny.*fqd + nx.*fpd;
end

return