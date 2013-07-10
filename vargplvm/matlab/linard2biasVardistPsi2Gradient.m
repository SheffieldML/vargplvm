function [gKern1, gKern2, gVarmeans, gVarcovars, gInd] = linard2biasVardistPsi2Gradient(linardKern, biasKern, vardist, Z, covGrad, learnInducing)

% LINARD2BIASVARDISTPSI2GRADIENT description.
  
% VARGPLVM
  
if nargin < 6
    learnInducing = 1;
end

% variational means
N = size(vardist.means,1);
%  inducing variables 
[M Q] = size(Z); 

[Psi2, Pnobias, Psi1] = linard2biasVardistPsi2Compute(linardKern, biasKern, vardist, Z);


% inverse variances
A = linardKern.inputScales;

% gradient for the bias parameter  
gKern2 = sum(sum(Pnobias.*covGrad)); 
Bnm = biasKern.variance*ones(size(Psi1)); 
BPsi1Covg = (Bnm*covGrad);

gInd = zeros(M, Q);

for q=1:vardist.latentDimension
   % 
   gKern(q) = sum(sum((vardist.means(:,q)*(Z(:,q)')).*BPsi1Covg));
   
   gVarmeans(:,q) = A(q)*sum((ones(vardist.numData,1)*Z(:,q)').*BPsi1Covg,2);
   if learnInducing
        gInd(:,q) = A(q)*sum((vardist.means(:,q)*ones(1,size(Z,1))).*BPsi1Covg,1)';
   end
   %   
end
%

gKern1 = 2*gKern(:)';  
% gVarmeans is N x Q matrix (N:number of data, Q:latent dimension)
% this will unfold this matrix column-wise 
gVarmeans = 2*gVarmeans(:)'; 

% gInd is M x Q matrix (M:number of inducing variables, Q:latent dimension)
% this will unfold this matrix column-wise 
gInd = 2*gInd(:)'; 

gVarcovars = zeros(1,prod(size(vardist.covars))); 


