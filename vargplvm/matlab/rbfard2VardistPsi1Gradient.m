function [gKern, gVarmeans, gVarcovars, gInd] = rbfard2VardistPsi1Gradient(rbfard2Kern, vardist, Z, covGrad, learnInducing)

% RBFARD2VARDISTPSI1GRADIENT description.
  
% VARGPLVM
  

if nargin < 5
    learnInducing = 1;
end


% variational means
N = size(vardist.means,1);
%  inducing variables 
M = size(Z,1); 


% evaluate the kernel matrix 
[K_fu Knovar] = rbfard2VardistPsi1Compute(rbfard2Kern, vardist, Z);

% inverse variances
A = rbfard2Kern.inputScales;

% gradient wrt variance of the kernel 
gKernvar = sum(sum(Knovar.*covGrad));  

KfuCovGrad = K_fu.*covGrad;

%--- new: preallocation
gVarmeans = zeros(vardist.numData, vardist.latentDimension);
gVarcovars = gVarmeans;
gInd = zeros(size(Z));
%--

% compute the gradient wrt lengthscales, variational means and variational variances  
for q=1:vardist.latentDimension
%
    S_q = vardist.covars(:,q);  
    Mu_q = vardist.means(:,q); 
    Z_q = Z(:,q)'; 
    
    % B3_q term (without A(q); see report)
    %B_q = repmat(1./(A(q)*S_q + 1), [1 M]).*(repmat(Mu_q,[1 M]) - repmat(Z_q,[N 1]));
    B_q = (repmat(Mu_q,[1 M]) - repmat(Z_q,[N 1]))./repmat(A(q)*S_q + 1, [1 M]);
    
    % derivatives wrt variational means and inducing inputs 
    %tmp = A(q)*((K_fu.*B_q).*covGrad);
    tmp = (B_q.*KfuCovGrad);
    
    % variational means: you sum out the columns (see report)
    gVarmeans(:,q) = -A(q)*sum(tmp,2); 
    
    if learnInducing
        % inducing inputs: you sum out the rows 
        gInd(:,q) = A(q)*sum(tmp,1)'; 
    end
    
    % 
    %B_q = repmat(1./(A(q)*S_q + 1), [1 M]).*dist2(Mu_q, repmat(Z_q);
    B_q = (B_q.*(repmat(Mu_q,[1 M]) - repmat(Z_q,[N 1])));
    
    % B1_q term (see report)
    %B1_q = -(0.5./repmat((A(q)*S_q + 1), [1 M])).*(repmat(S_q, [1 M]) + B_q);
    B1_q = (repmat(S_q, [1 M]) + B_q)./repmat((A(q)*S_q + 1), [1 M]);
    
    % gradients wrt kernel hyperparameters (lengthscales) 
    %gKernlengcs(q) = sum(sum((K_fu.*B1_q).*covGrad)); 
    gKernlengcs(q) = -0.5*sum(sum(B1_q.*KfuCovGrad)); 
    
    % B2_q term (see report)  
    %B1_q = ((0.5*A(q))./repmat((A(q)*S_q + 1), [1 M])).*(B_q - 1); 
  
    %B1_q = ((0.5*A(q))./repmat((A(q)*S_q + 1), [1 M])).*(A(q)*B_q - 1); 
    
    % gradient wrt variational covars (diagonal covariance matrices) 
    %gVarcovars(:,q) = sum((K_fu.*B1_q).*covGrad,2);
    gVarcovars(:,q) = sum((KfuCovGrad./repmat((A(q)*S_q + 1), [1 M])).*(A(q)*B_q - 1),2);
    
    %
end
%

gKern = [gKernvar gKernlengcs];

% gVarmeans is N x Q matrix (N:number of data, Q:latent dimension)
% this will unfold this matrix column-wise 
%gVarmeans = gVarmeans'; 
gVarmeans = gVarmeans(:)'; 

% gVarcovars is N x Q matrix (N:number of data, Q:latent dimension)
% this will unfold this matrix column-wise 
%gVarcovars = gVarcovars'; 
gVarcovars = 0.5*repmat(A,[N 1]).*gVarcovars;
gVarcovars = gVarcovars(:)';

% gInd is M x Q matrix (M:number of inducing variables, Q:latent dimension)
% this will unfold this matrix column-wise 
%gInd = gInd'; 
gInd = gInd(:)'; 

