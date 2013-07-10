function g = vargplvmPointGradient(x, model, y)
% VARGPLVMPOINTGRADIENT Wrapper function for gradient of a single point.
% FORMAT
% DESC is a wrapper function for the gradient of the log likelihood
% with respect to a point in the latent space. The GP-LVM
% model is one that is assumed to have already been trained.
% ARG x : the position in the latent space that is being optimised.
% ARG model : the trained GP-LVM model that is being optimised.
% ARG y : the position in data space for which the latent point is
% being optimised.
% RETURN g : the gradient of the log likelihood with respect to the
% latent position.
%
% SEEALSO : vargplvmPointLogLikeGradient, vargplvmOptimisePoint
%
% COPYRIGHT Michalis K. Titsias and Neil D. Lawrence, 2009, 2011
% COPYRIGHT Andreas C. Damianou, 2011

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
g = - vargplvmPointLogLikeGradient(model, vardistx, y);
