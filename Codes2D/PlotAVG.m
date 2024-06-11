function PlotAVG(x, y, u,S,b,interp,TRI)

% function [TRI,xout,yout,uout,interp] = PlotField2D(Nout, xin, yin, uin)
% Purpose: filled contour plot of solution data

% Globals2D;



uSave=u;
u(S.vmapM)=(u(S.vmapM)+u(S.vmapP))/2; 
u(S.ind)=repelem(accumarray(S.id,uSave(S.ind))./S.counts,S.counts);


PlotField2D(x, y, u,b,interp,TRI);
return
