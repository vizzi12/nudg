function rhsU=SWERHS2DViscosityTest(U,BC,BC2,fluxtype,t,S,S2,gauss,cub,viscositytype)
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
% dUdx= pagemtimes(cub.Dr',cub.W.*(cub.rx.*cU))+pagemtimes(cub.Ds',cub.W.*(cub.sx.*cU));
% dUdy= pagemtimes(cub.Dr',cub.W.*(cub.ry.*cU))+pagemtimes(cub.Ds',cub.W.*(cub.sy.*cU));

switch viscositytype
    case {'DB'}
        mu = ViscosityDB(U,S);
    case {'MDA'}
        mu = ViscosityMDA(U,gU,S,gauss);
    case {'EV'}
        [mu,S] = ViscosityEV(U,S,cub);
end

mu=ViscositySmooth(mu,S,S2);

dUdx= pagemtimes(S.V*S.V'*(cub.Dr)'.*cub.W(:,1)',(cub.rx.*cU))+pagemtimes(S.V*S.V'*(cub.Ds)'.*cub.W(:,1)',(cub.sx.*cU));
dUdy= pagemtimes(S.V*S.V'*(cub.Dr)'.*cub.W(:,1)',(cub.ry.*cU))+pagemtimes(S.V*S.V'*(cub.Ds)'.*cub.W(:,1)',(cub.sy.*cU));


fluxU=(gUM+gUP)/2;
% W1=-dUdx+pagemtimes(S.V*S.V'*gauss.interp'.*gauss.W(:,1)',gauss.nx.*gauss.sJ.*fluxU./gauss.J);
% W2=-dUdy+pagemtimes(S.V*S.V'*gauss.interp'.*gauss.W(:,1)',gauss.ny.*gauss.sJ.*fluxU./gauss.J);

drx=cub.J.*(cub.rx.*cU);
dry=cub.J.*(cub.ry.*cU);
dsx=cub.J.*(cub.sx.*cU);
dsy=cub.J.*(cub.sy.*cU);

surfacex=gauss.nx.*gauss.sJ.*fluxU;
surfacey=gauss.ny.*gauss.sJ.*fluxU;
rhsU=zeros(size(U));

W1=zeros(size(U));
W2=zeros(size(U));
for m=1:S.K
    k = m;
    drxk = pagemtimes((cub.mm(:,:,k)\cub.Dr').*cub.W(:,1)',drx(:,k,:));
    dryk = pagemtimes((cub.mm(:,:,k)\cub.Dr').*cub.W(:,1)',dry(:,k,:));
    dsxk = pagemtimes((cub.mm(:,:,k)\cub.Ds').*cub.W(:,1)',dsx(:,k,:));
    dsyk = pagemtimes((cub.mm(:,:,k)\cub.Ds').*cub.W(:,1)',dsy(:,k,:));
    surfacexk = pagemtimes((cub.mm(:,:,k)\gauss.interp').*gauss.W(:,1)', surfacex(:,k,:));
    surfaceyk = pagemtimes((cub.mm(:,:,k)\gauss.interp').*gauss.W(:,1)', surfacey(:,k,:));

    W1(:,k,:) = -drxk-dsxk+surfacexk;
    W2(:,k,:) = -dryk-dsyk+surfaceyk;
end

% W1(:,S.straight,:)=pagemtimes(S.V*S.V',(W1(:,S.straight,:)./S.J(:,S.straight)));
% W2(:,S.straight,:)=pagemtimes(S.V*S.V',(W2(:,S.straight,:)./S.J(:,S.straight)));

% Ncurved = length(S.curved);
% for m=1:Ncurved
%     k = S.curved(m);
%     mmCHOL = cub.mmCHOL(:,:,k);
%     W1(:,k,:) = pagemldivide(mmCHOL,pagemldivide(mmCHOL',W1(:,k,:)));
%     W2(:,k,:) = pagemldivide(mmCHOL,pagemldivide(mmCHOL',W2(:,k,:)));
% end


cW1=pagemtimes(cub.V,W1);
cW2=pagemtimes(cub.V,W2);

gW1 = pagemtimes(gauss.interp,W1);
gW2 = pagemtimes(gauss.interp,W2);

gW1M=gW1(gauss.mapM3D);
gW1P=gW1(gauss.mapP3D);

gW2M=gW2(gauss.mapM3D);
gW2P=gW2(gauss.mapP3D);

[gW1P,gW2P]=feval(BC2, gauss.x, gauss.y, gauss.nx, gauss.ny, gauss.mapW, gauss.mapD,gauss.mapO,gW1P,gW2P,t);

fluxW=gauss.nx.*(gW1M+gW1P)/2+gauss.ny.*(gW2M+gW2P)/2;

% ddrW = pagemtimes(S.V*S.V'*(cub.Dr)'.*cub.W(:,1)',(cub.rx.*cW1 + cub.ry.*cW2));
% ddsW = pagemtimes(S.V*S.V'*(cub.Ds)'.*cub.W(:,1)',(cub.sx.*cW1 + cub.sy.*cW2));
% 
% ddr = pagemtimes(S.V*S.V'*(cub.Dr)'.*cub.W(:,1)',(cub.rx.*F + cub.ry.*G));
% dds = pagemtimes(S.V*S.V'*(cub.Ds)'.*cub.W(:,1)',(cub.sx.*F + cub.sy.*G));
% rhsU = ddr + dds-ddrW -ddsW - pagemtimes(S.V*S.V'*gauss.interp'.*gauss.W(:,1)',gauss.sJ.*flux./gauss.J)+pagemtimes(S.V*S.V'*gauss.interp'.*gauss.W(:,1)',(gauss.sJ.*fluxW./gauss.J));

drW=cub.J.*(cub.rx.*cW1 + cub.ry.*cW2);
dsW=cub.J.*(cub.sx.*cW1 + cub.sy.*cW2);
dr=cub.J.*(cub.rx.*F + cub.ry.*G);
ds=cub.J.*(cub.sx.*F + cub.sy.*G);

surfaceW=gauss.sJ.*fluxW;
surface=gauss.sJ.*flux;

rhsU=zeros(size(U));
for m=1:S.K
    k = m;
    drWk = pagemtimes((cub.mm(:,:,k)\cub.Dr').*cub.W(:,1)',drW(:,k,:));
    dsWk = pagemtimes((cub.mm(:,:,k)\cub.Ds').*cub.W(:,1)',dsW(:,k,:));
    drk = pagemtimes((cub.mm(:,:,k)\cub.Dr').*cub.W(:,1)',dr(:,k,:));
    dsk = pagemtimes((cub.mm(:,:,k)\cub.Ds').*cub.W(:,1)',ds(:,k,:));
    surfaceWk = pagemtimes((cub.mm(:,:,k)\gauss.interp').*gauss.W(:,1)', surfaceW(:,k,:));
    surfacek = pagemtimes((cub.mm(:,:,k)\gauss.interp').*gauss.W(:,1)', surface(:,k,:));

    rhsU(:,k,:) = drk+dsk-dsWk-drWk-surfacek+surfaceWk;
end


% rhsU(:,S.straight,:)=pagemtimes(S.V*S.V',(rhsU(:,S.straight,:)./S.J(:,S.straight)));

% Ncurved = length(S.curved);
% for m=1:Ncurved
%     k = S.curved(m);
%     mmCHOL = cub.mmCHOL(:,:,k);
%     rhsU(:,k,:) = pagemldivide(mmCHOL,pagemldivide(mmCHOL',rhsU(:,k,:)));
% end

rhsU=rhsU+S.source(S.x,S.y,t,U);
return;


