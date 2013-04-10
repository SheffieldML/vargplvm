function [gKern, gVarmeans, gVarcovars, gInd] = rbfardVardistPsi2Gradient(rbfard2Kern, vardist, Z, covGrad)

% RBFARDVARDISTPSI2GRADIENT description.
  
% VARGPLVM
  
% variational means
N = size(vardist.means,1);
%  inducing variables 
M = size(Z,1); 


% evaluate the kernel matrix 
[K_fu Knovar] = rbfard2XvardistKernCompute(rbfard2Kern, vardist, Z);

% inverse variances
A = rbfard2Kern.inverseWidth;

% gradient wrt variance of the kernel 
gKernvar = sum(sum(Knovar.*covGrad));  

% compute the gradient wrt lengthscales, variational means and variational variances  
for q=1:vardist.latentDimension
%
    S_q = vardist.covars(:,q);  
    Mu_q = vardist.means(:,q); 
    Z_q = Z(:,q)'; 
    
    % B3_q term (see report)
    B_q = repmat(A(q)/(A(q)*S_q + 1), [1 M]).*(repmat(Mu_q,[1 M]) - repmat(Z_q,[N 1]));
 
    % derivatives wrt variational means and inducing inputs 
    tmp = (K_fu.*B_q).*covGrad;
    
    % variational means: you sum out the columns (see report)
    gVarmeans(:,q) = sum(tmp,2); 
    
    % inducing inputs: you sum out the rows 
    gInd(:,q) = sum(tmp,1)'; 
    
    %
    B_q = B_q.*(repmat(Mu_q,[1 M]) - repmat(Z_q,[N 1]));
    
    % B1_q term (see report)
    B1_q = -(0.5./repmat((A(q)*S_q + 1), [1 M])).*(repmat(S_q, [1 M]) + B_q);
    
    % gradients wrt kernel hyperparameters (lengthscales) 
    gKernlengcs(q) = sum(sum((K_fu.*B1_q).*covGrad)); 
    
    % B2_q term (see report) 
    B1_q = ((0.5*A(q))./repmat((A(q)*S_q + 1), [1 M])).*(B_q - 1); 
    
    % gradient wrt variational covars (diagonal covariance matrices) 
    gVarcovars(:,q) = sum((K_fu.*B1_q).*covGrad),2);
  
    %
end
     
gKern = [gKernlengcs gKernvar]; 