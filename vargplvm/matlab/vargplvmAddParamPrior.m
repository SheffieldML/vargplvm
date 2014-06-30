function model = vargplvmAddParamPrior(model, paramName, priorName, varargin)
% VARGPLVMADDPARAMPRIOR Add a prior in one of the parameters of the given
% vargplvm model
% DESC Add a prior in one of the parameters of the model.  This prior will
% NOT be optimised (ie fixed parameters).
% ARG model: The model for which we want to add a prior
% ARG paramName: The name of the parameter(s) for which the prior will be
% added (a regexp will match the given name)
% ARG priorName: The prior to be used. Make sure that the function
% [priorName 'priorParamInit'] exists.
% ARG varargin: Any arguments needed for the specific constructor of the
% prior are going to be pased
% 
% EXAMPLE: 
% % Add a prior on the noise parameter to encourage high SNR. High SNR is
% % achieved when beta is high compared to the data variance. So we will
% % constrain it relative to that value.
% >> varData = var(model.m(:));
% >> meanP = 400; stdP = 100;
% % Mean and std of the prior depend on the variance of the data, this says
% % that no matter what the variance of the data is, the SNR is encouraged to
% % be around meanP. 
% >> prec = 1/(stdP^2*varData);
% >> mu = meanP/varData;
% % prec is actually multiplied by the gradient. Should be selected so that
% % the overall gradient coming from this precision is comparable to the
% % gradient without the prior.
% >> model = vargplvmAddParamPrior(model, 'beta', 'gaussian2', [prec mu]);
% % Optional, for stronger effect
% >> model.paramPriors{1}.prior.scale = 100;
%
% % For example, for the gaussian2 prior, the gradient is:
% g = -(prior.precision*length(x)).*(x - prior.mean);
% so we select precision to be:
% g = vargplvmLogLikeGradients(model);
% ch = transformedParamChainRule(model.betaTransform, 1, model.beta);
% prec = g(model.paramPriors{1}.index)/(ch*model.paramPriors{1}.prior.dim*(model.beta-mu))
% model.paramPriors{1}.prior.precisions = abs(prec);
%
% SEEALSO: vargplvmParamPriorLogProb.m, vargplvmParamPriorGradients.m,
% priorParamInit.m
% 
% COPYRIGHT: Andreas C. Damianou, 2013
%
% VARGPLVM

[~, names] = vargplvmExtractParam(model);

if ~isfield(model, 'paramPriors')
    model.paramPriors = {};
end

%--- INDEX
switch paramName
    case {'beta', 'input scale'}
        tmpIndex = cellfun(@(x)regexpi(x,paramName),names,'un',0);
    otherwise
        error('Prior requested for unrecognised parameter')
end

paramPriors.index = [];
for i=1:length(tmpIndex)
    if ~isempty(tmpIndex{i})
        paramPriors.index = [paramPriors.index i];
    end
end
fprintf('# Added a %s prior on parameter: %s (%d indices)\n', priorName, paramName, length(paramPriors.index));

%---- PARAM. NAME
paramPriors.paramName = paramName;

%--- Actual prior
paramPriors.prior = priorCreate(priorName, varargin{:});

% Add prior
model.paramPriors{end+1} = paramPriors;

