function model = svargplvmExpandParam(model, params)

% SVARGPLVMEXPANDPARAM Expand a parameter vector into a shared variational GP-LVM model.
% FORMAT
% DESC takes a SVARGPLVM structure and a vector of parameters, and
% fills the structure with the given parameters. Also performs any
% necessary precomputation for likelihood and gradient
% computations, so can be computationally intensive to call.
% Parallelism can be applied w.r.t the models, so that some speed-up is achieved on
% multi-core machines.

% Parameters must be passed as a vector in the following order (left to right) 
% - parameter{size} -
% vardistParams{model.vardist.nParams} % mu, S
%       OR
% [dynamicsVardistParams{dynamics.vardist.nParams} dynamics.kernParams{dynamics.kern.nParams}] % mu_bar, lambda
%    % Followed by:
% private params of 1st model {model.comp{1}.nPrivateParams}
% private params of 2nd model {model.comp{2}.nPrivateParams}
%           ...
% private params of i-th model {model.comp{i}.nPrivateParams}
%
% ARG model : the svargplvm model to update with parameters
% ARG params : parameter vector
% RETURN model : model with updated parameters
%
% SEEALSO : svargplvmCreate, svargplvmExtractParam, modelExpandParam
%
% COPYRIGHT : Andreas C. Damianou, 2011

% VARGPLVM


try
    pool_open = matlabpool('size')>0;
catch e
    pool_open = 0;
end

if pool_open && (isfield(model,'parallel') && model.parallel)
    model = svargplvmExpandParamPar(model, params);
else
    model = svargplvmExpandParamOrig(model, params);
end
%!!! Since all dynamical models have the same "dynamics" structure, we can
%avoid the expandParam for ALL the dynamics structures. In fact, we should
%not store all these dynamics structures at all.


function model = svargplvmExpandParamPar(model, params)
%startVal = 1;
%endVal = model.N*model.q;
%paramsX = params(startVal:endVal);
% model.barmu = reshape(params(startVal:endVal), model.N, model.q);
startValS=1;

if isfield(model, 'dynamics') & ~isempty(model.dynamics)
    endValS = model.dynamics.nParams;
else
    endValS = model.vardist.nParams;
end
sharedParams = params(startValS:endValS);

startVal{1} = endValS+1;
endVal{1} = startVal{1} + model.comp{1}.nPrivateParams-1;
for i = 2:model.numModels
    % model.comp{i}.X = model.X;
    startVal{i} = endVal{i-1}+1;
    endVal{i} = startVal{i} + model.comp{i}.nPrivateParams-1; 
end

modelTempComp = model.comp;
parfor i = 1:model.numModels
    % model.comp{i}.X = model.X;
   % startVal = endVal+1;
   % endVal = startVal + model.comp{i}.nPrivateParams-1; 
    params_i = [sharedParams params(startVal{i}:endVal{i})];
    modelTempComp{i} = vargplvmExpandParam(modelTempComp{i}, params_i);
    if strcmp(modelTempComp{i}.kern.type, 'rbfardjit') || ~iscell(modelTempComp{i}.kern)
        modelTempComp{i}.kern.comp{1}.inputScales = modelTempComp{i}.kern.inputScales;
    end
    %[par, names] = vargplvmExtractParam(model.comp{i}); %%%%5
    %for j=5:-1:1 %%
    %params_i(end-j) %%%
    %end%%
    %model.comp{i}.kern.comp{1}.inputScales
end
model.comp = modelTempComp;


% with vargplvmExpandParam for the submodels, and given that the dynamics
% parameters (if there are dynamics) are shared, barmu has been turned into
% mu and this is the same for all submodels. The same holds for the
% variational distribution.
model.X = model.comp{1}.X;

if isfield(model, 'dynamics') & ~isempty(model.dynamics)
    model.dynamics = model.comp{1}.dynamics; % = model.comp{i}.dynamics for all submodels m
end

model.vardist = model.comp{1}.vardist;




function model = svargplvmExpandParamOrig(model, params)

%startVal = 1;
%endVal = model.N*model.q;
%paramsX = params(startVal:endVal);
% model.barmu = reshape(params(startVal:endVal), model.N, model.q);

startVal=1;

if isfield(model, 'dynamics') & ~isempty(model.dynamics)
    endVal = model.dynamics.nParams;
else
    endVal = model.vardist.nParams;
end
sharedParams = params(startVal:endVal);



for i = 1:model.numModels
    % model.comp{i}.X = model.X;
    startVal = endVal+1;
    endVal = startVal + model.comp{i}.nPrivateParams-1; 
    params_i = [sharedParams params(startVal:endVal)];
    model.comp{i} = vargplvmExpandParam(model.comp{i}, params_i);
    % For compatibility ...
    if strcmp(model.comp{i}.kern.type, 'rbfardjit')
        model.comp{i}.kern.comp{1}.inputScales = model.comp{i}.kern.inputScales;
    end
    %[par, names] = vargplvmExtractParam(model.comp{i}); %%%%5
    %for j=5:-1:1 %%
    %params_i(end-j) %%%
    %end%%
    %model.comp{i}.kern.comp{1}.inputScales
end

% with vargplvmExpandParam for the submodels, and given that the dynamics
% parameters (if there are dynamics) are shared, barmu has been turned into
% mu and this is the same for all submodels. The same holds for the
% variational distribution.
model.X = model.comp{1}.X;

if isfield(model, 'dynamics') & ~isempty(model.dynamics)
    model.dynamics = model.comp{1}.dynamics; % = model.comp{i}.dynamics for all submodels m
end

model.vardist = model.comp{1}.vardist;

