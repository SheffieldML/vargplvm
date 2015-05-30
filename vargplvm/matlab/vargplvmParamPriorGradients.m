function g = vargplvmParamPriorGradients(model, gs)
% VARGPLVMPARAMPRIORGRADIENTS Gradients of any priors put on parameters of
% a vargplvm model
% DESC Gradients of any priors put on parameters of
% a vargplvm model
% SEEALSO: vargplvmAddParamPrior.m, vargplvmParamPriorLogProb.m,
% priorParamInit.m
% COPYRIGHT: Andreas C. Damianou, 2013
% VARGPLVM

if nargin < 2
    if isfield(model, 'numParams')
        g = zeros(1, max(model.nParams, model.numParams));
    else
        % This should be the correct one
        g = zeros(1, model.nParams);
    end
else
    g = zeros(1, gs);
end
if ~isfield(model, 'paramPriors')
    return
else
    paramPriors = model.paramPriors;
end

for i=1:length(paramPriors)
    cur_params = paramPriors{i}.paramName;
    cur_ind = paramPriors{i}.index; % That's only for gradients
    switch cur_params
        case 'beta'
           set_param = model.beta;
           paramChainRule = transformedParamChainRule(model.betaTransform, 1, model.beta);
        case 'input scale'
            if ~isfield(model.kern, 'comp') || strcmp(model.kern.type, 'rbfardjit')
                set_param = model.kern.inputScales;
                paramChainRule = transformedParamChainRule(model.kern.transforms.type, model.kern.transforms.index, set_param);
            else
                set_param = [];
                paramChainRule = [];
                for cc = 1:length(model.kern.comp)
                    switch model.kern.comp{cc}.type
                        case {'rbfard','rbfard2','rbfardjit'}
                            set_param = [set_param model.kern.comp{cc}.inputScales];
                            paramChainRule = [paramChainRule transformedParamChainRule(model.kern.comp{cc}.transforms.type, model.kern.comp{cc}.transforms.index, set_param)];
                    end
                end
            end
        case 'dynamicsWhiteKernelVariance'
            set_param = [];
            paramChainRule = [];
            for cc = 1:length(model.dynamics.kern.comp)
                if strcmp(model.dynamics.kern.comp{cc}.type, 'white')
                    set_param = [set_param model.dynamics.kern.comp{cc}.variance];
                    paramChainRule = [paramChainRule transformedParamChainRule(model.dynamics.kern.comp{cc}.transforms.type, model.dynamics.kern.comp{cc}.transforms.index, set_param)];
                end
            end
    end
    cur_grad = priorGradient(paramPriors{i}.prior, set_param).*paramChainRule;
    if isfield(paramPriors{i}.prior, 'scale')
        cur_grad = cur_grad*paramPriors{i}.prior.scale;
    end
    g(cur_ind) = g(cur_ind) + cur_grad;
end