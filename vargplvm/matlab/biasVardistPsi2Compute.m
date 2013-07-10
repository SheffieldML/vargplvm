function [Psi2 P] = biasVardistPsi2Compute(biaskern, vardist, Z)

% BIASVARDISTPSI2COMPUTE one line description
% FORMAT
% DESC description
% RETURN Psi2 : description
% RETURN P : description
% ARG biasKern : the kernel structure associated with the bias kernel.
% ARG vardist : description
% ARG Z : description
%
%
% SEEALSO : others
%
%
% COPYRIGHT : Michalis K. Titsias, 2009
%

% VARGPLVM


Psi2 = repmat(vardist.numData*(biaskern.variance^2),size(Z,1),size(Z,1)); 
P = [];




