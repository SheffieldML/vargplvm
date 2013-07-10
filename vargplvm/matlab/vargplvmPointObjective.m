function f = vargplvmPointObjective(x, model, y)

% VARGPLVMPOINTOBJECTIVE Wrapper function for objective of a single point in latent space and the output location..
% FORMAT
% DESC provides a wrapper function for the negative log probability
% of a given data point under the posterior distribution of the
% Gaussian process induced by the training data.
% ARG x : location in input space for the point.
% ARG model : the model structure for which the negative log
% probability of the given data under the posterior is to be computed.
% ARG y : the location in data space for the point.
% RETURN f : the negative of the log probability of the given data
% point under the posterior distribution induced by the training data.
%
% SEEALSO : vargplvmCreate, vargplvmPointLogLikelihood, vargplvmOptimisePoint
%
% COPYRIGHT : Michalis K. Titsias and Neil D. Lawrence, 2009

% VARGPLVM

if isfield(model, 'dynamics') && ~isempty(model.dynamics)
    if isfield(model.dynamics, 'reoptimise') && model.dynamics.reoptimise  %%% RE-OPT-CODE-NEW
        [vardistx, model] = vargplvmPartExpand(model, x); %%% RE-OPT-CODE-NEW
    else %%% RE-OPT-CODE-NEW
        % this is doing the expand
        x = reshape(x, model.N+size(y,1), model.dynamics.q*2);
        xtrain = x(1:model.N,:);
        xtest = x(model.N+1:end,:);
        model.dynamics.vardist = vardistExpandParam(model.dynamics.vardist, xtrain);
        vardistx = vardistExpandParam(model.vardistx, xtest);
        % end of expand
    end %%% RE-OPT-CODE-NEW
else
    vardistx = model.vardistx;
    vardistx = vardistExpandParam(vardistx, x);
end
f = - vargplvmPointLogLikelihood(model, vardistx, y);