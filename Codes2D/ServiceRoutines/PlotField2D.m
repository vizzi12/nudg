function PlotField2D(xout, yout, uin,b,interp,TRI)

% function [TRI,xout,yout,uout,interp] = PlotField2D(Nout, xin, yin, uin)
% Purpose: filled contour plot of solution data

% Globals2D;
     
% build equally spaced grid on reference triangle


% interpolate node coordinates and field to equally spaced nodes
 uout = interp*uin;

% render and format solution field
if b
    h=trisurf(TRI, xout(:), yout(:), uout(:),'FaceColor',[0.3010 0.7450 0.9330],'faceAlpha',1,'Edgecolor','[0 0 0]','EdgeAlpha',0);
%     h.LineWidth=0.000000000000001;
%     set(h,'edgelighting','gouraud')
    lighting gouraud 
    camlight headlight,
else
    trisurf(TRI, xout(:), yout(:), uout(:),'FaceColor',[0.6250 0.3203 0.1758],'Edgecolor','None');
%     set(h,'edgelighting','gouraud')
end
% shading interp
% material shiny,    

return
