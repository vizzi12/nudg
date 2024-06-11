function [usol] = Advec2Dsol(x, y, t, lam)
    usol = exp(- (1/(2*lam^2)) * ( (x - 0.5*cos(pi-2*pi*t)).^2 + (y - 0.5*sin(pi-2*pi*t)).^2 ) );
end