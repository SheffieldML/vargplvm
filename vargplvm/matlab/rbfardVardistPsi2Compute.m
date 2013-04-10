function [K, Knovar, argExp] = rbfardVardistPsi2Compute(rbfardKern, vardist, Z)

% RBFARDVARDISTPSI2COMPUTE one line description
% FORMAT
% DESC description
% RETURN K : description
% RETURN Knovar : description
% RETURN argExp : description
% ARG rbfardKern : the kernel structure associated with the rbfard kernel.
% ARG vardist : description
% ARG Z : description
%
%
% SEEALSO : others
%
%
% COPYRIGHT : Michalis K. Titsias, 2009
%

% VARGPLVM

% variational means
N  = size(vardist.means,1);
%  inducing variables 
M = size(Z,1); 

A = rbfardKern.inverseWidth;
         
argExp = zeros(N,M); 
normfactor = ones(N,1);
for q=1:vardist.latentDimension
%
    S_q = vardist.covars(:,q);  
    normfactor = normfactor.*(A(q)*S_q + 1);
    Mu_q = vardist.means(:,q); 
    Z_q = Z(:,q)'; 
    distan = (repmat(Mu_q,[1 M]) - repmat(Z_q,[N 1])).^2;
    argExp = argExp + repmat(A(q)/(A(q)*S_q + 1), [1 M]).*distan;
%
end
Knovar = repmat(normfactor,[1 M]).*exp(-0.5*argExp); 
K = rbfardKern.variance.*Knovar; 


