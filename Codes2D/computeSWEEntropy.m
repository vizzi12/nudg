function S = computeSWEEntropy(U,S)
g=9.81;

h=U(:,:,1);
u=U(:,:,2)./h;
v=U(:,:,3)./h;

S.E=0.5*g*h.^2+0.5*u.^2.*h+0.5*v.^2.*h;
S.F(:,:,1)=g*u.*h.^2 + 0.5*u.^3.*h + 0.5*u.*v.^2.*h;
S.F(:,:,2)=g*v.*h.^2 + 0.5.*u.^2.*v.*h + 0.5*v.^3.*h;

dFdr= pagemtimes(S.Dr,S.F);
dFds= pagemtimes(S.Ds,S.F);

S.div=S.rx.*dFdr(:,:,1) + S.sx.*dFds(:,:,1)+S.ry.*dFdr(:,:,2) + S.sy.*dFds(:,:,2);
end