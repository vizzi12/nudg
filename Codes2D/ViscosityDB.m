function mu = ViscosityDB(U,S)
g=9.81;
c = sqrt(g*U(:,:,1));
V=U(:,:,2:3)./U(:,:,1);
dVdr= pagemtimes(S.Dr,V);
dVds= pagemtimes(S.Ds,V);
cBeta=0.5;
muBeta=cBeta*abs((S.rx.*dVdr(:,:,1) + S.sx.*dVds(:,:,1)+(S.ry.*dVdr(:,:,2) + S.sy.*dVds(:,:,2)))).*(S.dtscale/S.N).^2;
cMax=1;
muMax=cMax*(S.dtscale/S.N).*max(sqrt(sum(V.^2,3))+c);





mu=min(muBeta,muMax.*ones(S.Np,1));

