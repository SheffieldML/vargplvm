function f = varDimSequenceObjective(params,model,X,dim,varargin)

% VARDIMSEQUENCEOBJECTIVE Objective of output variance over
% subspace sequence
% FORMAT
% DESC Computes outputvariance associated with latent sequence
% ARG params : latent subspace sequence
% ARG model : fgplvm model generating observation
% ARG X : full latent sequence
% ARG dim : latent dimensions to alter
% RETURN f : objective
%
% SEEALSO : varDimSequenceGradient
%
% COPYRIGHT : Neil D. Lawrence, Carl Henrik Ek, 2008

% SGPLVM

X(:,dim) = reshape(params,size(X,1),length(dim));
Y = gpPosteriorMeanVar(model,X);
X = X(:)';

f = fgplvmSequenceObjective(X,model,Y,varargin{:});

return;
