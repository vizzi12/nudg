function [mu,S] = ViscosityEV(U,S,cub)
g=9.81;
c = sqrt(g*U(:,:,1));

h=U(:,:,1);
V=U(:,:,2:3)./h;


R=(S.E-S.Eold)/S.dt+0.5*(S.div+S.divold);
if anynan(R)
    R=zeros(size(R));
end

FM=S.F(S.mapM3D);
FP=S.F(S.mapP3D);

Fjump=(S.N./S.dtscale).*(S.nx.*(FM(:,:,1)-FP(:,:,1))+S.ny.*(FM(:,:,2)-FP(:,:,2)));

int=sum(cub.W.*(cub.V*S.E),'all');

T=max(abs(S.E(:)-1/S.omega*int));
% if any(R > 0)
%     2+2;
% end

D=max(max(abs(R)),max(abs(Fjump)));

D=D./T;

cE=0.15;
mu=cE*(S.dtscale/S.N).^2.*D;


cMax=1;
muMax=cMax*(S.dtscale/S.N).*max(sqrt(sum(V.^2,3))+c);
mu=min(mu,muMax.*ones(S.Np,1));