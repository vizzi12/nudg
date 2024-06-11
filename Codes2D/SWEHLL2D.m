function flux = SWEHLL2D(nx, ny, UM, UP,form)
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

FLh= qMr;
FRh= qPr;

fh=( sR.*FLh-sL.*FRh+sR.*sL.*(hP-hM) ) ./ (sR - sL);
fh(sL >= 0) = FLh(sL >= 0);
fh(sR <= 0) = FRh(sR <= 0);

FLq= qMr.*uMr+0.5*g*hM.^2;
FRq= qPr.*uPr+0.5*g*hP.^2;

fq=( sR.*FLq-sL.*FRq+sR.*sL.*(qPr-qMr) ) ./ (sR - sL);
fq(sL >= 0) = FLq(sL >= 0);
fq(sR <= 0) = FRq(sR <= 0);

FLp= uMr.*vMr.*hM;
FRp= uPr.*vPr.*hP;

fp=( sR.*FLp-sL.*FRp+sR.*sL.*(pPr-pMr) ) ./ (sR - sL);
fp(sL >= 0) = FLp(sL >= 0);
fp(sR <= 0) = FRp(sR <= 0);

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