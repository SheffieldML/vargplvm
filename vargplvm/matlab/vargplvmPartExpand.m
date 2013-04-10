function [vardistx, model] = vargplvmPartExpand(model, x, update)
% VARGPLVMPARTEXPAND
% This function is used when the model is using jointly the test and
% training data (possibly with missing values from the first) to optimise
% and in that case only part of the parameters are being optimised. This
% function, thus, expands the model by taking into account only the
% optimised parameters (the old ones are being kept fixed).

% In this version of the function the params that are considered to be
% optimised again are the var. distribution, dyn. kernel hyperparams and
% inducing points. The rest (base kernel and beta) are being kept fixed.

% x = [mu_bar lambda theta_t X_u] where mu_bar and lambda are augmented to
% include the vardist. of the test points.

% Note on the order of the parameters:
% [(N*Q*2)              |theta_t| Q*k |theta_f| 1] -> vargplvmExtractParam(model)
% [(N*Q*2)+(Nstar*Q*2)  |theta_t| Q*k]             -> x
%
% COPYRIGHT: Andreas C. Damianou, Michalis K. Titsias 2011
% VARGPLVM

% The model params before optimisation
paramsOrig = vargplvmExtractParam(model);

% The optimised parameters x here are not only the vardist. but also
% gDynKern and inducing points.
paramsNew = x;

% x will now hold only the vardist. params.
vardistxNparams=length(model.dynamics.t_star) * model.q * 2;
x = paramsNew(1:(model.dynamics.vardist.nParams + vardistxNparams));

% now separate the variational disribution into the training part and the
% testing part and update the original training model (only with the new training
% variational distribution) and the test variational distribution
% this is doing the expand
vardistNumData = size(x,2)/(model.q*2);
x = reshape(x, vardistNumData, model.dynamics.q*2);
xtrain = x(1:model.N,:);
xtest = x(model.N+1:end,:);
%model.dynamics.vardist = vardistExpandParam(model.dynamics.vardist, xtrain);
vardistx = vardistExpandParam(model.vardistx, xtest(:)');

% The model will be expanded with the origParams vector where all the
% new values (optimised) will overwrite the old ones.
params = paramsOrig;

startVal = 1;
endVal = model.dynamics.vardist.nParams;

% This will overwrite the vardist original params with the newly optimised
% ones
params(startVal:endVal) = xtrain(:)';

% This will overwrite the dyn.kern params and the inducing points.
startVal = endVal+1;
endVal = endVal + model.dynamics.kern.nParams + model.q * model.k;
params(startVal:endVal) = paramsNew((startVal+vardistxNparams):(endVal+vardistxNparams));

% Update the model
if exist('update') && update
    model = vargplvmExpandParam(model, params);
else
    % Update the model withoug calling update stats (this is used in the
    % optimisation likelihood and wrapper functions which do all the
    % precomputations themselves)
    model = vargplvmExpandParamNoUpdate(model, params);
end


%-----------


% This is copy-paste of vargplvmExpandParam without the last call to
% vargplvmUpdateStats (because this is done inside the objective and gradient functions)
function model = vargplvmExpandParamNoUpdate(model,params)
startVal = 1;
if isfield(model, 'dynamics') & ~isempty(model.dynamics)
    % variational parameters (reparametrized) AND dyn.kernel's parameters
    endVal = model.dynamics.nParams;
    model.dynamics = modelExpandParam(model.dynamics, params(startVal:endVal));
else
    % variational parameters (means and covariances), original ones
    endVal = model.vardist.nParams;
    model.vardist = modelExpandParam(model.vardist, params(startVal:endVal));
end

% inducing inputs
startVal = endVal+1;
if model.fixInducing
    if isfield(model, 'dynamics') & ~isempty(model.dynamics)
        Kt = kernCompute(model.dynamics.kern, model.dynamics.t);%%%%%%%%%%%%5
        model.X_u = Kt*model.dynamics.vardist.means; %dynamics
        model.X_u = model.X_u(model.inducingIndices,:);
    else
        model.X_u = model.vardist.means(model.inducingIndices, :); % static
    end
    % X_u values are taken from X values.
    % model.X_u = model.X(model.inducingIndices, :);
else
    % Parameters include inducing variables.
    endVal = endVal + model.q*model.k;
    model.X_u = reshape(params(startVal:endVal),model.k,model.q);
end

% kernel hyperparameters
startVal = endVal+1;
endVal = endVal + model.kern.nParams;
model.kern = kernExpandParam(model.kern, params(startVal:endVal));

% likelihood beta parameters
if model.optimiseBeta
    startVal = endVal + 1;
    endVal = endVal + prod(size(model.beta));
    fhandle = str2func([model.betaTransform 'Transform']);
    model.beta = fhandle(params(startVal:endVal), 'atox');
end

model.nParams = endVal;
