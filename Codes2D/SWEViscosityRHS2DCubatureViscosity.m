function rhsU=SWEViscosityRHS2DCubatureViscosity(U,BC,BC2,fluxtype,t,S,S2,gauss,cub,viscositytype)
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

[F,G] = SWEFluxes2D(cU);
dUdx= pagemtimes(cub.Dr',cub.W.*(cub.rx.*cU))+pagemtimes(cub.Ds',cub.W.*(cub.sx.*cU));
dUdy= pagemtimes(cub.Dr',cub.W.*(cub.ry.*cU))+pagemtimes(cub.Ds',cub.W.*(cub.sy.*cU));

fluxU=(gUM+gUP)/2;
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
  case {'EB'}
    mu = ViscosityEB(U,S);
  case {'MDA'}
    mu = ViscosityMDA(U,gU,S,gauss);
end

mu=ViscositySmooth(mu,S,S2);


W1=mu.*W1;
W2=mu.*W2;

cu=cU(:,:,2)./cU(:,:,1);
cv=cU(:,:,3)./cU(:,:,1);
dudx= (cub.Dr')*(cub.W.*(cub.rx.*cu))+(cub.Ds')*(cub.W.*(cub.sx.*cu));
dvdy= (cub.Dr')*(cub.W.*(cub.ry.*cv))+(cub.Ds')*(cub.W.*(cub.sy.*cv));

nu=0.01252836781;
W1(:,:,2)=W1(:,:,2)+nu*U(:,:,1).*dudx;
W2(:,:,3)=W2(:,:,3)+nu*U(:,:,1).*dvdy;


cW1=pagemtimes(cub.V,W1);
cW2=pagemtimes(cub.V,W2);

gW1 = pagemtimes(gauss.interp,W1);
gW2 = pagemtimes(gauss.interp,W2);

gW1M=gW1(gauss.mapM3D);
gW1P=gW1(gauss.mapP3D);

gW2M=gW2(gauss.mapM3D);
gW2P=gW2(gauss.mapP3D);

[gW1P,gW2P]=feval(BC2, gauss.x, gauss.y, gauss.nx, gauss.ny, gauss.mapW, gauss.mapD,gauss.mapO,gW1P,gW2P);

fluxW=gauss.nx.*(gW1M+gW1P)/2+gauss.ny.*(gW2M+gW2P)/2;

ddrW = pagemtimes(cub.Dr',cub.W.*(cub.rx.*cW1 + cub.ry.*cW2));
ddsW = pagemtimes(cub.Ds',cub.W.*(cub.sx.*cW1 + cub.sy.*cW2));

ddr = pagemtimes(cub.Dr',cub.W.*(cub.rx.*F + cub.ry.*G));
dds = pagemtimes(cub.Ds',cub.W.*(cub.sx.*F + cub.sy.*G));
rhsU = ddr + dds-ddrW -ddsW - pagemtimes(gauss.interp',(gauss.W.*flux))+pagemtimes(gauss.interp',(gauss.W.*fluxW));
rhsU(:,S.straight,:)=pagemtimes(S.V*S.V',(rhsU(:,S.straight,:)./S.J(:,S.straight)));

Ncurved = length(S.curved);
for m=1:Ncurved
    k = S.curved(m);
    mmCHOL = cub.mmCHOL(:,:,k);
    rhsU(:,k,:) = pagemldivide(mmCHOL,pagemldivide(mmCHOL',rhsU(:,k,:)));
end

rhsU=rhsU+S.source(S.x,S.y,t,U);
return;


