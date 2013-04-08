function g = varDimGradient(params,model,X,dim)

% VARDIMGRADIENT Output variance gradient
% FORMAT
% DESC Compute latent gradients for output space variance
% ARG params : latent locations for optimised dimensions
% ARG model : fgplvm model generating observation space
% ARG X : complete latent location
% ARG dim : optimised dimensions
% RETURN g : gradients
%
% SEEALSO : varDimObjective
%
% COPYRIGHT : Neil D. Lawrence, Carl Henrik Ek, 2008
%
% Copied from the SGPLVM toolbox.
% 
% SHEFFIELDML

X(:,dim) = params;

[void g] = gpPosteriorGradMeanVar(model,X);
g = g(dim,1)';

return;