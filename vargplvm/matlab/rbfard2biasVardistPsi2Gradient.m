function [gKern1, gKern2, gVarmeans, gVarcovars, gInd] = rbfard2biasVardistPsi2Gradient(rbfardKern, biasKern, vardist, Z, covGrad, learnInducing)

% RBFARD2BIASVARDISTPSI2GRADIENT description.
  
% VARGPLVM
  
if nargin < 6
    learnInducing = 1;
end

% variational means
N = size(vardist.means,1);
%  inducing variables 
[M Q] = size(Z); 

[Psi2, Pnobias, Psi1] = rbfard2biasVardistPsi2Compute(rbfardKern, biasKern, vardist, Z);


% inverse variances
A = rbfardKern.inputScales;

% gradient wrt variance of the rbfard2 kernel 
gKernvar = sum(sum(Psi2.*covGrad))/rbfardKern.variance;  
% gradient for the bias parameter  
gKern2 = sum(sum(Pnobias.*covGrad)); 
Bnm = biasKern.variance*ones(size(Psi1)); 
BPsi1Covg = Psi1.*(Bnm*covGrad); 

%-- New: preallocation
    gVarmeans = zeros(N,Q);
    gVarcovars = zeros(N,Q);
    gInd = zeros(M,Q);
%---


% compute the gradient wrt lengthscales, variational means and variational variances  
for q=1:vardist.latentDimension
%
    S_q = vardist.covars(:,q);  
    Mu_q = vardist.means(:,q); 
    Z_q = Z(:,q)'; 
    
    % B3_q term (without A(q); see report)
    B_q = (repmat(Mu_q,[1 M]) - repmat(Z_q,[N 1]))./repmat(A(q)*S_q + 1, [1 M]);
    
    % derivatives wrt variational means and inducing inputs 
    tmp = (B_q.*BPsi1Covg);
    
    % variational means: you sum out the columns (see report)
    gVarmeans(:,q) = -A(q)*sum(tmp,2); 
    
    % inducing inputs: you sum out the rows 
    if learnInducing
        gInd(:,q) = A(q)*sum(tmp,1)'; 
    end
    
    % 
    B_q = (B_q.*(repmat(Mu_q,[1 M]) - repmat(Z_q,[N 1])));
    
    % B1_q term (see report)
    B1_q = (repmat(S_q, [1 M]) + B_q)./repmat((A(q)*S_q + 1), [1 M]);
    
    % gradients wrt kernel hyperparameters (lengthscales) 
    gKernlengcs(q) = -0.5*sum(sum(B1_q.*BPsi1Covg)); 
    
    % gradient wrt variational covars (diagonal covariance matrices) 
    gVarcovars(:,q) = sum((BPsi1Covg./repmat((A(q)*S_q + 1), [1 M])).*(A(q)*B_q - 1),2);
    
    %
end
%

gKern1 = [gKernvar 2*gKernlengcs];

% gVarmeans is N x Q matrix (N:number of data, Q:latent dimension)
% this will unfold this matrix column-wise
gVarmeans = 2*gVarmeans(:)'; 

% gVarcovars is N x Q matrix (N:number of data, Q:latent dimension)
% this will unfold this matrix column-wise 
gVarcovars = repmat(A,[N 1]).*gVarcovars;
gVarcovars = gVarcovars(:)';

% gInd is M x Q matrix (M:number of inducing variables, Q:latent dimension)
% this will unfold this matrix column-wise 
gInd = 2*gInd(:)'; 


