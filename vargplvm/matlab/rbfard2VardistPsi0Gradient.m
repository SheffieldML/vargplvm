function [gKern, gVarmeans, gVarcovars] = rbfard2VardistPsi0Gradient(rbfard2Kern, vardist, covGrad)

% RBFARD2VARDISTPSI0GRADIENT Description
  
% VARGPLVM
gKern = zeros(1,rbfard2Kern.nParams); 
gKern(1) = covGrad*vardist.numData;
 
gVarmeans = zeros(1,prod(size(vardist.means))); 
gVarcovars = zeros(1,prod(size(vardist.means))); 




