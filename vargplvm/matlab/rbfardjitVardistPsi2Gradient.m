function [gKern, gVarmeans, gVarcovars, gInd] = rbfardjitVardistPsi2Gradient(rbfardKern, vardist, Z, covGrad, learnInducing)

% RBFARDJITVARDISTPSI2GRADIENT description.
  
% VARGPLVM
  
if nargin < 5
    learnInducing = 1;
end

[gKern, gVarmeans, gVarcovars, gInd] = rbfard2VardistPsi2Gradient(rbfardKern, vardist, Z, covGrad, learnInducing);
