function f = varDimObjective(params,model,X,dim)

% VARDIMOBJECTIVE Objective of output variance over subspace
% FORMAT
% DESC Computes outputvariance associated with latent location
% ARG params : latent subspace
% ARG model : fgplvm model generating observation
% ARG X : full latent location
% ARG dim : latent dimensions to alter
% RETURN f : objective
%
% SEEALSO : varDimGradient
%
% COPYRIGHT : Neil D. Lawrence, Carl Henrik Ek, 2008
%
% Copied from the SGPLVM toolbox...
%
% SHEFFIELDML

X(:,dim) = params;

[void f] = gpPosteriorMeanVar(model,X);
f = f(1,1);

return;