function [U] = SWELimiter2Dn(U,BC,S,limit,t)
% function [LQ] = EulerLimiter2D(Q, SolutionBC, time)
% Purpose: limit the Euler solution using slope limiting adapted from 
% A SLOPE LIMITING PROCEDURE IN DISCONTINUOUS GALERKIN FINITE ELEMENT METHOD FOR 
% GASDYNAMICS APPLICATIONS. SHUANGZHANG TU AND SHAHROUZ ALIABADI 
% INTERNATIONAL JOURNAL OF NUMERICAL ANALYSIS AND MODELING, Volume 2, Number 2, Pages 163
% Globals2D;

% Gas constant
g = 9.81;


% 1. compute geometric information for 4 element patch containing each element
% Build average matrix

% 2. Find cell averages of conserved & primitive variables in each 4 element patch

% Compute cell averages of conserved variables
% hC = limit.AVE*U(:,:,1); qC = limit.AVE*U(:,:,2); pC = limit.AVE*U(:,:,3);
% 
% % Find neighbor values of conserved variables
% PC(:,:,1)=hC(EToE'); PC(:,:,2)=qC(EToE'); PC(:,:,3)=pC(EToE');

UC=pagemtimes(limit.AVE,U);
PC=UC(cat(3,S.EToE',S.EToE'+S.K,S.EToE'+2*S.K));

% Apply boundary conditions to cell averages of ghost cells at boundary faces
[PC] = feval(BC,[limit.xc1;limit.xc2;limit.xc3], [limit.yc1;limit.yc2;limit.yc3], limit.fnx, limit.fny,limit.idW, limit.idD, 0, PC,t);


Uavg=squeeze(UC)';

Ul=pagemtimes(S.invV,U); Ul(3:S.N+1,:,:)=0; Ul(S.N+3:S.Np,:,:)=0;

Vnew=Vandermonde2D(1, [0,0,-1]', [-1,0,0]');

Um=pagemtimes(Vnew,Ul([1,2,S.N+2],:,:));
Um=permute(Um,[3,2,1]);

Am=Um-Uavg;
Am=reshape(Am,[3,1,S.K,3]);

F=permute(PC,[3,2,1]);

Bm=permute(limit.alpha(:,1,:),[2,3,1 ]).*(F-Uavg)+permute(limit.alpha(:,2,:),[2,3,1]).*(F(:,:,[2,3,1])-Uavg);
Bm=reshape(Bm,[3,1,S.K,3]);

c=reshape(sqrt(g*UC(:,:,1)),1,1,[]);
u=reshape(UC(:,:,2)./UC(:,:,1),1,1,[]);
v=reshape(UC(:,:,3)./UC(:,:,1),1,1,[]);

sx=limit.s(:,1,:,:);
sy=limit.s(:,2,:,:);

L=1/(2*c) .* [u.*sx+v.*sy+c,-sx,-sy ; 2*c.*(-u.*sy+v.*sx),2*c.*sy,-2*c.*sx ; -u.*sx-v.*sy+c,sx,sy];
R=[ones(size(sx)),sx.*0,ones(size(sx)) ; u-c.*sx,sy,u+c.*sx ; v-c.*sy,-sx,v+c.*sy];

LA=pagemtimes(L,Am);
LB=1.5*pagemtimes(L,Bm);

delta=zeros(size(LA));
M=0.01;

delta(:)=minmodB([LA(:)';LB(:)'],M,S.hgrid);

delta=pagemtimes(R,delta);
U00=sum(delta,4);

idx= sum(abs(U00)) > 1e-10;
if any(idx)
   pos=sum(max(0,delta(:,:,idx,:)),4);
   neg=sum(max(0,-delta(:,:,idx,:)),4);
   thetaP=min(1,neg./pos); thetaM=min(1,pos./neg);
   delta(:,:,idx,:)=thetaP.*max(0,delta(:,:,idx,:))-thetaM.*max(0,-delta(:,:,idx,:));
end

delta=squeeze(delta)+Uavg;

idx = any(any(abs(delta-Um) > 1e-10,3));

if any(idx)
    Um(:,idx,:)=delta(:,idx,:);
    Um=permute(Um,[3,2,1]);
    
    Ul([1,2,S.N+2],idx,:)=pagemldivide(Vnew,Um(:,idx,:));
    U(:,idx,:)=pagemtimes(S.V,Ul(:,idx,:));
end

return;