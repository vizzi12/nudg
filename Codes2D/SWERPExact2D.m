function flux = SWERPExact2D(nx, ny, UM, UP,form)
g=9.81;

hM = UM(:,:,1); qM = UM(:,:,2); pM = UM(:,:,3);
hP = UP(:,:,1); qP = UP(:,:,2); pP = UP(:,:,3);

qMr = nx.*qM + ny.*pM;
pMr = -ny.*qM + nx.*pM;

qPr = nx.*qP + ny.*pP;
pPr = -ny.*qP + nx.*pP;

uMr=qMr./hM;
vMr=pMr./hM;

uPr=qPr./hP;
vPr=pPr./hP;

fh=zeros(size(uMr));
fq=zeros(size(uMr));
fp=zeros(size(uMr));
for i=1:length(uMr(:))
%     if (hM(i)-hP(i) > 4.9)
%         2+2;
%     end
    [hLOCAL,uLOCAL,vLOCAL]=SWRPExact(hM(i),uMr(i),vMr(i),hP(i),uPr(i),vPr(i),1,0,0,1);
    fh(i)=hLOCAL*uLOCAL;
%     if uLOCAL > 2
%         2+2;
%     end
    fq(i)=hLOCAL*uLOCAL^2+0.5*g*hLOCAL^2;
    fp(i)=hLOCAL*uLOCAL*vLOCAL;
end


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