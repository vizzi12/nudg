function [F] = Filter1D(NN,Nc,s)

% function [F] = Filter1D(N,Nc,s)
% Purpose : Initialize 1D filter matrix of size N.
%           Order of exponential filter is (even) s with cutoff at Nc;

Globals1D;
filterdiag = ones(NN+1,1);
alpha = -log(eps);

% Initialize filter function
for i=Nc:NN
    filterdiag(i+1) = exp(-alpha*((i-Nc)/(NN-Nc))^s);
end;

F = V*diag(filterdiag)*invV;
return;
