function k0 = rbfardjitVardistPsi0Compute(rbfardKern, vardist)

% RBFARDJITVARDISTPSI0COMPUTE description.
  
% SHEFFIELDML
  
% variational means


% dont intlude the jitter term (it only for the inducing matrix K_uu)

k0 = vardist.numData*rbfardKern.variance; 
