function mu = ViscositySmooth(mu,S,S2)

Vnew=Vandermonde2D(S.N,S2.r,S2.s);
muModal=S.invV*mu;
mu=Vnew*muModal;

muSave=mu;
mu(S2.vmapM)=(mu(S2.vmapM)+mu(S2.vmapP))/2; 

% mu(1,:)=(muSave(1,:)+muSave(S2.vmapP(1,:))+muSave(S2.vmapP(7,:)))/3;
% mu(3,:)=(muSave(1,:)+muSave(S2.vmapP(3,:))+muSave(S2.vmapP(4,:)))/3;
% mu(6,:)=(muSave(1,:)+muSave(S2.vmapP(9,:))+muSave(S2.vmapP(6,:)))/3;

mu(S2.ind)=repelem(accumarray(S2.id,muSave(S2.ind))./S2.counts,S2.counts);

muModal=S2.invV*mu;

Vnew2=Vandermonde2D(S2.N,S.r,S.s);
mu=Vnew2*muModal;
return