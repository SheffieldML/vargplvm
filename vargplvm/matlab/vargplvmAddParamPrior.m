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

