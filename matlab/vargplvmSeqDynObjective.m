function f = vargplvmSeqDynObjective(x, model, y)

% VARGPLVMSEQDYNOBJECTIVE Wrapper function for objective of a group of points in latent space and the output locations..
% FORMAT
% DESC provides a wrapper function for the negative log probability
% of a group of data points under the posterior distribution of the
% Gaussian process induced by the training data.
% ARG x : locations in input space for the point.
% ARG model : the model structure for which the negative log
% probability of the given data under the posterior is to be computed.
% ARG y : the location in data space for the points.
% RETURN f : the negative of the log probability of the given data
% point under the posterior distribution induced by the training data.
% 
% SEEALSO : vargplvmCreate, vargplvmPointLogLikelihood, vargplvmOptimisePoint
%
% COPYRIGHT : Michalis K. Titsias and Neil D. Lawrence, 2011

% VARGPLVM

% % this is doing the expand 
% x = reshape(x, model.N+size(y,1), model.dynamics.q*2);
% xtrain = x(1:model.N,:);
% xtest = x(model.N+1:end,:);
% model.dynamics.vardist = vardistExpandParam(model.dynamics.vardist, xtrain);
% vardistx = vardistExpandParam(model.vardistx, xtest);
% % end of expand 
vardistx = model.vardistx;
vardistx = vardistExpandParam(vardistx, x);

f = - vargplvmSeqDynLogLikelihood(model, vardistx, y);
%f = -TEMPvargplvmPointLogLikelihoodSTATIC(model, vardistx, y);