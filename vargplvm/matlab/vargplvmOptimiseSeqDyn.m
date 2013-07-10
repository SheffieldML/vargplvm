function [X varX] = vargplvmOptimiseSeqDyn(model, vardistx, y, display, iters);

% VARGPLVMOPTIMISESEQDYN Optimise the positions of a group of latent
% point which correspond to a independent sequence.
% FORMAT
% DESC optimises the locations of a group of points in latent space
% given an initialisation and a group of (partly) observed data points. 
% ARG model : the model for which the points will be optimised.
% ARG vardistx : the initialisation of the points in the latent space.
% ARG y : the observed data point for which the latent points are to
% be optimised.
% ARG display : whether or not to display the iterations of the
% optimisation (default: true)
% ARG iters : maximum number of iterations for the optimisation
% (default 2000).
% RETURN x : the optimised means in the latent space.
% RETURN varx : the optimised variances in the latent space.

%
% COPYRIGHT :  Michalis K. Titsias 2011
% SEEALSO : vargplvmCreate, vargplvmOptimiseSequence, vargplvmPointObjective, vargplvmPointGradient

% VARGPLVM

if nargin < 5
  iters = 2000;
  %if nargin < 5
    display = true;
  %end
end

options = optOptions;
if display
  options(1) = 1;
  %options(9) = 1;
end
options(14) = iters;


if isfield(model, 'optimiser')
  optim = str2func(model.optimiser);
else
  optim = str2func('scg');
end


% % augment the training and testing variational distributions
% vardist = vardistCreate(zeros(model.N+size(y,1), model.q), model.q, 'gaussian');
% vardist.means = [model.dynamics.vardist.means; vardistx.means];
% vardist.covars = [model.dynamics.vardist.covars; vardistx.covars]; 
% vardist.numData = size(vardist.means,1);
% vardist.nParams = 2*prod(size(vardist.means));
x = modelExtractParam(vardistx);

if strcmp(func2str(optim), 'optimiMinimize')
  % Carl Rasmussen's minimize function 
  x = optim('vargplvmSeqDynObjectiveGradient', x, options, model, y);
else
  % NETLAB style optimization.
  x = optim('vargplvmSeqDynObjective', x,  options, ...
            'vargplvmSeqDynGradient', model, y);
end

% now separate the variational disribution into the training part and the
% testing part and update the original training model (only with the new training 
% variational distribution) and the test variational distribution
% this is doing the expand 

% x = reshape(x, vardist.numData, model.dynamics.q*2);
% xtrain = x(1:model.N,:);
% xtest = x(model.N+1:end,:);
% model.dynamics.vardist = vardistExpandParam(model.dynamics.vardist, xtrain);
% vardistx = vardistExpandParam(model.vardistx, xtest);
% % end of expand 

vardistx = vardistExpandParam(vardistx,x);
X = vardistx.means;
varX = vardistx.covars;
