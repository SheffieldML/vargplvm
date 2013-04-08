function [K, P] = biasVardistPsi1Compute(biaskern, vardist, Z)

% BIASVARDISTPSI1COMPUTE one line description
% FORMAT
% DESC description
% RETURN K : description
% RETURN P : description
% ARG biasKern : the kernel structure associated with the white kernel.
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

K = repmat(biaskern.variance,size(vardist.means,1),size(Z,1));

P = [];