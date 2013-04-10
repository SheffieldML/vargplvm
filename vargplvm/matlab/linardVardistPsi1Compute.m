function [K, Knovar, argExp] = linardVardistPsi1Compute(linardkern, vardist, Z)

% LINARDVARDISTPSI1COMPUTE description.

% VARGPLVM
  
% variational means
N  = size(vardist.means,1);
%  inducing variables 
M = size(Z,1); 

A = rbfardKern.inverseWidth;
  
scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;
    
if nargin < 3
  sk = x*x';
else
  x2 = x2*scales;
  sk = x*x2';
end
k = sk*kern.variance;       



