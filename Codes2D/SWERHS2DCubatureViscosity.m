function [rhsU]=SWERHS2DCubatureViscosity(U,BC,BC2,fluxtype,t,S,S2,gauss,cub,viscositytype)
g=9.81;

gU = pagemtimes(gauss.interp,U);

gUM=gU(gauss.mapM3D);
gUP=gU(gauss.mapP3D);

gUP=feval(BC, gauss.x, gauss.y, gauss.nx, gauss.ny, gauss.mapW, gauss.mapD,gauss.mapO,gUP,t);

switch fluxtype
  case {'HLL'}
    flux = SWEHLL2D(gauss.nx, gauss.ny, gUM, gUP,'weak'); 
  case {'HLLC'}
    flux = SWEHLLC2D(gauss.nx, gauss.ny, gUM, gUP,'weak'); 
  case {'RPExact'}
    flux = SWERPExact2D(gauss.nx, gauss.ny, gUM, gUP,'weak');
end

cU = pagemtimes(cub.V,U);

cUp = cU; 
cUp(:,:,2:3)=cUp(:,:,2:3)./cUp(:,:,1);

[F,G] = SWEFluxes2D(cU);
dUdx= pagemtimes(cub.Dr',cub.W.*(cub.rx.*cUp))+pagemtimes(cub.Ds',cub.W.*(cub.sx.*cUp));
dUdy= pagemtimes(cub.Dr',cub.W.*(cub.ry.*cUp))+pagemtimes(cub.Ds',cub.W.*(cub.sy.*cUp));


gUMp=gUM;
gUMp(:,:,2:3)=gUMp(:,:,2:3)./gUMp(:,:,1);

gUPp=gUP;
gUPp(:,:,2:3)=gUPp(:,:,2:3)./gUPp(:,:,1);
fluxU=(gUMp+gUPp)/2;
W1=-dUdx+pagemtimes(gauss.interp',(gauss.W.*gauss.nx.*fluxU));
W2=-dUdy+pagemtimes(gauss.interp',(gauss.W.*gauss.ny.*fluxU));

W1(:,S.straight,:)=pagemtimes(S.V*S.V',(W1(:,S.straight,:)./S.J(:,S.straight)));
W2(:,S.straight,:)=pagemtimes(S.V*S.V',(W2(:,S.straight,:)./S.J(:,S.straight)));

Ncurved = length(S.curved);
for m=1:Ncurved
    k = S.curved(m);
    mmCHOL = cub.mmCHOL(:,:,k);
    W1(:,k,:) = pagemldivide(mmCHOL,pagemldivide(mmCHOL',W1(:,k,:)));
    W2(:,k,:) = pagemldivide(mmCHOL,pagemldivide(mmCHOL',W2(:,k,:)));
end

switch viscositytype
    case {'DB'}
        mu = ViscosityDB(U,S);
    case {'MDA'}
        mu = ViscosityMDA(U,gU,S,gauss);
    case {'EV'}
        [mu,S] = ViscosityEV(U,S,cub);
end

mu=ViscositySmooth(mu,S,S2);


W1=mu.*W1;
W2=mu.*W2;


cW1=pagemtimes(cub.V,W1);
cW2=pagemtimes(cub.V,W2);

cW1(:,:,2:3)=cW1(:,:,2:3).*cU(:,:,1);
cW2(:,:,2:3)=cW2(:,:,2:3).*cU(:,:,1);

gW1 = pagemtimes(gauss.interp,W1);
gW2 = pagemtimes(gauss.interp,W2);


gW1M=gW1(gauss.mapM3D);
gW1P=gW1(gauss.mapP3D);

gW2M=gW2(gauss.mapM3D);
gW2P=gW2(gauss.mapP3D);

[gW1P,gW2P]=feval(BC2, gauss.x, gauss.y, gauss.nx, gauss.ny, gauss.mapW, gauss.mapD,gauss.mapO,gW1P,gW2P,t);

gW1M(:,:,2:3)=gW1M(:,:,2:3).*gUM(:,:,1);
gW2M(:,:,2:3)=gW2M(:,:,2:3).*gUM(:,:,1);

gW1P(:,:,2:3)=gW1P(:,:,2:3).*gUP(:,:,1);
gW2P(:,:,2:3)=gW2P(:,:,2:3).*gUP(:,:,1);

fluxW=gauss.nx.*(gW1M+gW1P)/2+gauss.ny.*(gW2M+gW2P)/2;

ddrW = pagemtimes(cub.Dr',cub.W.*(cub.rx.*cW1 + cub.ry.*cW2));
ddsW = pagemtimes(cub.Ds',cub.W.*(cub.sx.*cW1 + cub.sy.*cW2));

% ddrW(:,:,1)=0*ddrW(:,:,1);
% ddsW(:,:,1)=0*ddrW(:,:,1);
% 
% fluxW(:,:,1)=0*fluxW(:,:,1);

ddr = pagemtimes(cub.Dr',cub.W.*(cub.rx.*F + cub.ry.*G));
dds = pagemtimes(cub.Ds',cub.W.*(cub.sx.*F + cub.sy.*G));
rhsU = ddr + dds-ddrW -ddsW - pagemtimes(gauss.interp',(gauss.W.*flux))+pagemtimes(gauss.interp',(gauss.W.*fluxW));
% rhsU = -pagemtimes(gauss.interp',(gauss.W.*flux))+pagemtimes(gauss.interp',(gauss.W.*fluxW));
% rhsU = ddr + dds-ddrW -ddsW;

rhsU(:,S.straight,:)=pagemtimes(S.V*S.V',(rhsU(:,S.straight,:)./S.J(:,S.straight)));

Ncurved = length(S.curved);
for m=1:Ncurved
    k = S.curved(m);
    mmCHOL = cub.mmCHOL(:,:,k);
    rhsU(:,k,:) = pagemldivide(mmCHOL,pagemldivide(mmCHOL',rhsU(:,k,:)));
end

rhsU=rhsU+S.source(S.x,S.y,t,U);
return;


