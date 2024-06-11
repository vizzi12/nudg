function [ux,uy] = Grad2D(u, G);

% function [ux,uy] = Grad2D(u);
% Purpose: Compute 2D gradient field of scalar u

ur = G.Dr*u; us = G.Ds*u;
ux = G.rx.*ur + G.sx.*us; uy = G.ry.*ur + G.sy.*us;
return
