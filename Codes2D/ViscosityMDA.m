function mu = ViscosityMDA(U,gU,S,gauss)
g=9.81;
h=U(:,:,1);
V=U(:,:,2:3)./h;

c = sqrt(g*U(:,:,1));
cMax=1;
muMax=cMax*(S.dtscale/S.N).*max(sqrt(sum(V.^2,3))+c);

hhat(:,:,1)=S.V1D1\h(S.Fmask(:,1),:);
hhat(:,:,2)=S.V1D2\h(S.Fmask(:,2),:);
hhat(:,:,3)=S.V1D3\h(S.Fmask(:,3),:);

hhat=hhat(2:S.Nfp,:,:);

m=(1:S.N)';
b=m.^(-S.N)./sqrt(sum(m.^(-2*S.N)));

gh=gU(:,:,1);
% L2=gauss.sJ.*gh.^2.*gauss.W;
L2=gh.^2.*gauss.W;
hL2 = zeros(1,S.K,3);
hL2(1,:,1) = sum(L2(1:gauss.NGauss,:),1);
hL2(1,:,2) = sum(L2(gauss.NGauss+1:2*gauss.NGauss,:),1);
hL2(1,:,3) = sum(L2(2*gauss.NGauss+1:3*gauss.NGauss,:),1);



hhat=sqrt(  hhat.^2+ pagemtimes(b.^2,hL2));

hhat(S.N,:,:) = max(hhat(S.N-1:S.N,:,:));
hhat = cummax(hhat,1,"reverse");

A=[ones(S.N,1),-log(m)];
tau=pagemldivide(A,log(hhat));
tau=tau(2,:,:);

tau=min(tau,[],3);

mu=tau*0+(tau<1).*muMax+muMax.*(1 <= tau & tau <= 3).*(1-(tau-1)/2);
mu=ones(S.Np,1)*mu;