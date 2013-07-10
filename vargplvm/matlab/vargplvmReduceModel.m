function [mm, order] = vargplvmReduceModel(model, P, dims)

% VARGPLVMREDUCEMODEL prunes out dimensions of the model.
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
% MODIFICATIONS : Neil D. Lawrence, 2009, Patrick Sauer, 2011
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

mm.vardist.means  = model.vardist.means(:,order(1:P));
mm.vardist.covars = model.vardist.covars(:,order(1:P));

mm.X_u = model.X_u(:,order(1:P));
mm.X   = model.vardist.means(:,order(1:P));
mm.inputSclsOrder = order;

initParams   = vargplvmExtractParam(mm);
mm.numParams = length(initParams);

% This forces kernel computation.
mm = vargplvmExpandParam(mm, initParams);
