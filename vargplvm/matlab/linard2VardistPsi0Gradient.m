function [gKern, gVarmeans, gVarcovars] = linard2VardistPsi0Gradient(linard2Kern, vardist, covGrad)

% LINARD2VARDISTPSI0GRADIENT description.

% VARGPLVM
  
A = linard2Kern.inputScales;
gKern = covGrad*sum((vardist.means.*vardist.means) + vardist.covars,1); 
 
gVarmeans = 2*(vardist.means*sparse(diag(A))); 
%gVarmeans1 = 2*(repmat(A,size(vardist.means,1),1).*vardist.means); 

gVarcovars = ones(size(vardist.means,1),1)*A; 

gVarmeans = covGrad*gVarmeans(:)'; 
gVarcovars = covGrad*gVarcovars(:)';


