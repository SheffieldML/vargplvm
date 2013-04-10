function [f, g] = vargplvmPointObjectiveGradient(x, model, y)

% VARGPLVMPOINTOBJECTIVEGRADIENT Wrapper function for objective and gradient of a single point in latent space and the output location..
% FORMAT
% DESC provides a wrapper function for the negative log probability
% of a given data point under the posterior distribution of the
% Gaussian process induced by the training data. Also returns the
% gradient of the negative log probability with respect to the
% given latent point.
% ARG x : location in input space for the point.
% ARG model : the model structure for which the negative log
% probability of the given data under the posterior is to be computed.
% ARG y : the location in data space for the point.
% RETURN f : the negative of the log probability of the given data
% point under the posterior distribution induced by the training data.
% RETURN g : the gradient of the log probability with respect to
% the given latent point.
% 
% SEEALSO : vargplvmCreate, vargplvmPointLogLikelihood,
% vargplvmOptimisePoint, vargplvmObjective, vargplvmGradient
%
% COPYRIGHT : Michalis K. Titsias and Neil D. Lawrence, 2009

% VARGPLVM

% Check how the optimiser has given the parameters
if size(xvec, 1) > size(xvec, 2)
  % As a column vector ... transpose everything.
  transpose = true;
  x = x';
else
  transpose = false;
end
f = - vargplvmPointLogLikelihood(model, x, y);

if nargout > 1
  g = - vargplvmPointLogLikeGradient(model, x, y);
end
if transpose
  g = g';
end

