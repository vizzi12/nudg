function [IM] = InterpMatrix2D(rout, sout,S)

% function [IM] = InterpMatrix2D(rout, sout)
% Purpose: Compute local elemental interpolation matrix

% Globals2D;
 
% compute Vandermonde at (rout,sout)
Vout = Vandermonde2D(S.N, rout, sout);

% build interpolation matrix
IM = Vout*S.invV;
return

  
  
