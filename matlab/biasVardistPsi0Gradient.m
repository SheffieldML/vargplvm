function [gKern, gVarmeans, gVarcovars] = biasVardistPsi0Gradient(biaskern, vardist, covGrad)

% BIASVARDISTPSI0GRADIENT one line description
% FORMAT
% DESC description here.
% RETURN gKern :
% RETURN gVarmeans :
% RETURN gVarcovars :
% RETURN gInd :
% ARG biaskern : the kernel structure associated with the bias kernel.
% ARG vardist :
% ARG Z : 
% ARG covGrad : 
%
%
% SEEALSO : others
%
%
% COPYRIGHT : Michalis K. Titsias, 2009
%

% VARGPLVM

gKern = covGrad*vardist.numData;
 
gVarmeans = zeros(1,prod(size(vardist.means))); 
gVarcovars = zeros(1,prod(size(vardist.means))); 