function g = varDimSequenceGradient(params,model,X,dim,varargin)

% VARDIMSEQUENCEGRADIENT Output variance gradient for latent sequence
% FORMAT
% DESC Compute latent gradients for output space variance of a
% sequence
% ARG params : latent locations for optimised dimensions
% ARG model : fgplvm model generating observation space
% ARG X : complete latent location
% ARG dim : optimised dimensions
% RETURN g : gradients
%
% SEEALSO : varDimSequenceObjective
%
% COPYRIGHT : Neil D. Lawrence, Carl Henrik Ek, 2008
%
% Copied from the SGPLVM toolbox...
%
% VARGPLVM

X(:,dim) = reshape(params,size(X,1),length(dim));
Y = gpPosteriorMeanVar(model,X);
X = X(:)';

g = fgplvmSequenceGradient(X,model,Y,varargin{:});
g = reshape(g,size(Y,1),model.q);
g = g(:,dim);
g = g(:)';

return;