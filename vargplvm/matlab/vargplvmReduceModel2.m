function [mm, order] = vargplvmReduceModel2(model, P, dims)

% VARGPLVMREDUCEMODEL2 prunes out dimensions of the model.
% FORMAT
% DESC order the latent dimensions acorrding to the inputScales and
% reduces the model to have smaller number of latent dimensions.
% ARG model : the model to be reduced.
% ARG P : the number of dimensions to move to (setting to model.q will
% just reorder the dimensions in the model).
% ARG dims : (optional) explicit set of dimensions to use
% RETURN model : the model with the reduced number of dimensions.
%
% COPYRIGHT : Michalis K. Titsias, 2009
%  
% MODIFICATIONS : Neil D. Lawrence, 2009, Patrick Sauer, 2011, Andreas
% Damianou 2012
% 
% SEEALSO : vargplvmCreate  

% VARGPLVM

if nargin == 3
    P = length(dims);
end


% create temporary model
options = vargplvmOptions('dtcvar');
options.kern =[];
options.numActive = model.k;
if ~strcmp(model.kern.type,'cmpnd')
  options.kern = model.kern.type;
elseif strcmp(model.kern.type,'cmpnd')
  options.kern{1} = model.kern.comp{1}.type;
  for i = 2:length(model.kern.comp)
    options.kern{i} = model.kern.comp{i}.type;
  end
end

if isstruct( model.betaTransform )
    options.betaTransform = model.betaTransform;
end   

mm = vargplvmCreate(P, model.d, model.y, options);
N = size(model.vardist.means,1);

if ~strcmp(model.kern.type,'cmpnd')
  % 
  if strcmp(model.kern.type,'rbfardjit') | strcmp(model.kern.type,'linard2') | strcmp(model.kern.type,'rbfard2')
    %
    if nargin == 2
        [vals, order] = sort(-model.kern.inputScales);
        mm.kern.inputScales = model.kern.inputScales(order(1:P));
    else
        order = [dims setdiff(1:length(model.kern.inputScales), dims)];
        mm.kern.inputScales = model.kern.inputScales(order);
    end        
    %
  end
  %
else
  %
  for i = 1:length(model.kern.comp)
    %
    if strcmp(model.kern.comp{i}.type,'rbfardjit') | strcmp(model.kern.comp{i}.type,'linard2') | strcmp(model.kern.comp{i}.type,'rbfard2')
      %  
      if nargin == 2
          [vals, order] = sort(-model.kern.comp{i}.inputScales);
          mm.kern.comp{i}.inputScales = model.kern.comp{i}.inputScales(order(1:P));
      else
          order = [dims setdiff(1:length(model.kern.comp{i}.inputScales), dims)];
          mm.kern.comp{i}.inputScales = model.kern.comp{i}.inputScales(order);
      end
      % you order only wrt the first ARD kernel you find 
      break;  
      %
    end
    %
  end
  %
end


%-----
% This code segment is added so that the old parameters (the ones that
% still remain in the reduced model) are copied in the new model. With the
% old code (vargplvmReduceModel), the optimised models were "reset" to the
% default values for the parameters, apart from the lengthscales which were
% copied correctly. A test to see it works correcty is to check that:
% vargplvmLogLikelihood(model) - vargplvmLogLikelihood(vargplvmReduceModel2(model, model.q))
% is very small.

mm.vardist.means  = model.vardist.means(:,order(1:P));
mm.vardist.covars = model.vardist.covars(:,order(1:P));

mm.X_u = model.X_u(:,order(1:P));
mm.X   = model.vardist.means(:,order(1:P));
mm.inputSclsOrder = order;

%dimsToKeep = order(1:P);
%dimsToDiscard = setdiff(order,dimsToKeep);

% Vardist
mm.vardist.means  = model.vardist.means(:,order(1:P));
mm.vardist.covars = model.vardist.covars(:,order(1:P));
mm.vardist.latentDimension = P;
mm.vardist.nParams = mm.vardist.numData*mm.vardist.latentDimension*2;
if isfield(model.vardist, 'transforms') && length(model.vardist.transforms.index)==length(model.N*model.q+1:model.N*model.q*2)
    mm.vardist.transforms.index = mm.N * mm.q+1 : mm.N*mm.q*2;
end


% Dynamics...
if isfield(model, 'dynamics') && ~isempty(model.dynamics)
    mm.dynamics = model.dynamics;
    mm.dynamics.q = P;
    mm.dynamics.vardist.means  = model.dynamics.vardist.means(:,order(1:P));
    mm.dynamics.vardist.covars = model.dynamics.vardist.covars(:,order(1:P));
    mm.dynamics.vardist.latentDimension = P;
    mm.dynamics.vardist.nParams = mm.dynamics.vardist.numData*mm.dynamics.vardist.latentDimension*2;
    if isfield(model.dynamics.vardist, 'transforms') && length(model.dynamics.vardist.transforms.index)==length(model.N*model.q+1:model.N*model.q*2)
        mm.dynamics.vardist.transforms.index = mm.N * mm.q+1 : mm.N*mm.q*2;
    end
    mm.dynamics.nParams = model.dynamics.nParams - (model.dynamics.vardist.nParams - mm.dynamics.vardist.nParams);
end

% Parameters
params = [];

if isfield(mm, 'dynamics') & ~isempty(mm.dynamics)
    % [VariationalParameters(reparam)   dynKernelParameters]
    dynParams = modelExtractParam(mm.dynamics);
    params = dynParams;
else
    % Variational parameters 
    varParams = modelExtractParam(mm.vardist);
    params = varParams;
end

% Inducing inputs 
if ~model.fixInducing
    params =  [params mm.X_u(:)'];
end


% Kernel parameters  
[kernParams, names] = kernExtractParam(model.kern);
for i=1:length(names)
    if strfind(names{i}, 'input scale')
        scalesStart = i;
        break
    end
end
kernParamsNew=[];
scalesEnd = scalesStart + model.q-1;
scalesInd = scalesStart:scalesStart+model.q-1;
%scalesRem = setdiff(1:length(kernParams), scalesInd);
kernParamsNew = kernParams(1:scalesStart-1);
kernParamsNew = [kernParamsNew kernParams(scalesInd(order(1:P))) kernParams(scalesEnd+1:end)];
params = [params kernParamsNew];


% beta in the likelihood 
if model.optimiseBeta
   
   if ~isstruct(model.betaTransform)
       fhandle = str2func([model.betaTransform 'Transform']);
       betaParam = fhandle(model.beta, 'xtoa');
   else
      if isfield(model.betaTransform,'transformsettings') && ~isempty(model.betaTransform.transformsettings)
          fhandle = str2func([model.betaTransform.type 'Transform']);
          betaParam = fhandle(model.beta, 'xtoa', model.betaTransform.transformsettings);
      else
          error('vargplvmExtractParam: Invalid transform specified for beta.'); 
      end
   end   

   params = [params betaParam(:)'];
   
end
%-------


mm.numParams = length(params);

% This forces kernel computation.
mm = vargplvmExpandParam(mm, params);
