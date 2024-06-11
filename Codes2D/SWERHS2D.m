function [rhsU] = SWERHS2D(U,BC,fluxtype,form,t,S)
L=length(U(:))/3;
UM=U(cat(3,S.vmapM,S.vmapM+L,S.vmapM+2*L));
UP=U(cat(3,S.vmapP,S.vmapP+L,S.vmapP+2*L));

UP=feval(BC, S.Fx, S.Fy, S.nx, S.ny, S.mapW, S.mapD,S.mapO,UP,t);


switch fluxtype
  case {'HLL'}
    flux = SWEHLL2D(S.nx, S.ny, UM, UP,form); 
  case {'HLLC'}
    flux = SWEHLLC2D(S.nx, S.ny, UM, UP,form); 
end

[F,G] = SWEFluxes2D(U);

switch form
    case {'weak'}
        dFdr= pagemtimes(S.Drw,F);
        dGdr= pagemtimes(S.Drw,G);
        dFds= pagemtimes(S.Dsw,F);
        dGds= pagemtimes(S.Dsw,G);
        rhsU = (S.rx.*dFdr + S.sx.*dFds) + (S.ry.*dGdr + S.sy.*dGds) - pagemtimes(S.LIFT,(S.Fscale.*flux))+S.source(S.x,S.y,t,U);
    case {'strong'}
        dFdr= pagemtimes(S.Dr,F);
        dGdr= pagemtimes(S.Dr,G);
        dFds= pagemtimes(S.Ds,F);
        dGds= pagemtimes(S.Ds,G);
        rhsU = -(S.rx.*dFdr + S.sx.*dFds) -(S.ry.*dGdr + S.sy.*dGds) +pagemtimes(S.LIFT,(S.Fscale.*flux))+S.source(S.x,S.y,t,U);
end
return;


