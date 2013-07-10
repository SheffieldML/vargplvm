function [K, P] = whiteVardistPsi1Compute(whitekern, vardist, Z)

% WHITEVARDISTPSI1COMPUTE one line description
% FORMAT
% DESC description
% RETURN K : description
% RETURN P : description
% ARG whiteKern : the kernel structure associated with the white kernel.
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



  K = zeros(size(vardist.means,1), size(Z,1));
  
  P = [];