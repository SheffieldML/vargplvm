function [gKern, gVarmeans, gVarcovars, gInd] = rbfardjitVardistPsi1Gradient(rbfard2Kern, vardist, Z, covGrad, learnInducing)

% RBFARDJITVARDISTPSI1GRADIENT description.
  
% VARGPLVM
  
if nargin < 5
    learnInducing = 1;
end

 [gKern, gVarmeans, gVarcovars, gInd] = rbfard2VardistPsi1Gradient(rbfard2Kern, vardist, Z, covGrad, learnInducing);
