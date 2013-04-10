function [gKern, gVarmeans, gVarcovars, gInd] = biasVardistPsi1Gradient(biaskern, vardist, Z, covGrad, learnInducing)

% BIASVARDISTPSI1GRADIENT Compute gradient of bias variational PSI1.
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
% SEEALSO : 
%
% COPYRIGHT : Michalis K. Titsias, 2009
%

% VARGPLVM
  
gKern = sum(sum(ones(vardist.numData,size(Z,1)).*covGrad)); 

gVarmeans = zeros(1,prod(size(vardist.means))); 

gInd = zeros(1,prod(size(Z))); 

gVarcovars = zeros(1,prod(size(vardist.covars))); 
