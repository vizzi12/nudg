function rhsU=SWERHS2DCubature(U,BC,fluxtype,t,S,gauss,cub)
 
% gUM = zeros(gauss.NGauss*Nfaces, K, 3); gUP = zeros(gauss.NGauss*Nfaces, K, 3);
% for n=1:3
%   gU = gauss.interp*U(:,:,n);
%   gUM(:,:,n) = gU(gauss.mapM);  gUP(:,:,n) = gU(gauss.mapP);
% end 

gU = pagemtimes(gauss.interp,U);

gUM=gU(cat(3,gauss.mapM,gauss.mapM+length(gauss.mapM(:)),gauss.mapM+2*length(gauss.mapM(:))));
gUP=gU(cat(3,gauss.mapP,gauss.mapP+length(gauss.mapP(:)),gauss.mapP+2*length(gauss.mapP(:))));

gUP=feval(BC, gauss.x, gauss.y, gauss.nx, gauss.ny, gauss.mapW, gauss.mapD,gauss.mapO,gUP,t);

switch fluxtype
  case {'HLL'}
    flux = SWEHLL2D(gauss.nx, gauss.ny, gUM, gUP,'weak'); 
  case {'HLLC'}
    flux = SWEHLLC2D(gauss.nx, gauss.ny, gUM, gUP,'weak'); 
end

cU = pagemtimes(cub.V,U);

[F,G] = SWEFluxes2D(cU);

rhsU = zeros(S.Np, S.K, 3);
for n=1:3
  ddr = (cub.Dr')*(cub.W.*(cub.rx.*F(:,:,n) + cub.ry.*G(:,:,n)));
  dds = (cub.Ds')*(cub.W.*(cub.sx.*F(:,:,n) + cub.sy.*G(:,:,n)));
  rhsU(:,:,n) = ddr + dds - gauss.interp'*(gauss.W.*flux(:,:,n));
  rhsU(:,:,n)=S.V*S.V'*(rhsU(:,:,n)./S.J);
end

% ddr = pagemtimes(S.V*S.V'*(cub.Dr)'.*cub.W(:,1)',(cub.rx.*F + cub.ry.*G));
% dds = pagemtimes(S.V*S.V'*(cub.Ds)'.*cub.W(:,1)',(cub.sx.*F + cub.sy.*G));
% rhsU = ddr + dds - pagemtimes(S.V*S.V'*gauss.interp',(gauss.W.*flux)./gauss.J); %W=W*SJ
% rhsU(:,S.straight,:)=pagemtimes(,(rhsU(:,S.straight,:))); %./S.J(:,S.straight);

% ddr=pagemtimes((cub.Dr)'.*cub.W(:,1)',cub.J.*(cub.rx.*F + cub.ry.*G));
% dds=pagemtimes((cub.Ds)'.*cub.W(:,1)',cub.J.*(cub.sx.*F + cub.sy.*G));
% int=pagemtimes(gauss.interp',(gauss.W.*flux));
% rhsU = ddr + dds - pagemtimes(gauss.interp',(gauss.W.*flux));
% rhsU(:,S.curved,:) = ddr(:,S.curved,:) + dds(:,S.curved,:) - int(:,S.curved,:);

% ddr=cub.J.*(cub.rx.*F + cub.ry.*G);
% dds=cub.J.*(cub.sx.*F + cub.sy.*G);
% surface=gauss.sJ.*flux;
% rhsU=zeros(size(U));

% Ncurved = length(S.curved);
% for m=1:S.K
%     k = m;
%     ddrk = pagemtimes((cub.mm(:,:,k)\cub.Dr').*cub.W(:,1)',ddr(:,k,:));
%     ddsk = pagemtimes((cub.mm(:,:,k)\cub.Ds').*cub.W(:,1)',dds(:,k,:));
%     surfacek = pagemtimes((cub.mm(:,:,k)\gauss.interp').*gauss.W(:,1)', surface(:,k,:));
% 
%     rhsU(:,k,:) = ddrk+ddsk-surfacek; %
% end

% for m=1:S.K
%     k = m;
%     mmCHOL = cub.mmCHOL(:,:,k);
%     rhsU(:,k,:) = pagemldivide(mmCHOL,pagemldivide(mmCHOL',rhsU(:,k,:))); %
% %     rhsU(:,k,:) = pagemldivide(cub.mm(:,:k),rhsU(:,k,:));
% end

rhsU=rhsU+S.source(S.x,S.y,t,U);
return;


