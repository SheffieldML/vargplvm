function [params, names] = vardistExtractParam(vardist)

% VARDISTEXTRACTPARAM Extract a parameter vector from a vardist structure.
% FORMAT
% DESC extracts a parameter vector from a given VARDIST structure.
% ARG model : the model from which parameters are to be extracted.
% RETURN params : the parameter vector extracted from the model.
%
% DESC does the same as above, but also returns parameter names.
% ARG model : the model structure containing the information about
% the model.
% RETURN params : a vector of parameters from the model.
% RETURN names : cell array of parameter names.
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% SEEALSO : vardistCreate, vardistExpandParam, modelExtractParam

% VARGPLVM



%means = vardist.means'; 
%covs = vardist.covars';

% the variational means and diagonal covariances obtained COLUMN-WISE 
params = [vardist.means(:)' vardist.covars(:)'];

% names
if nargout > 1  
    for i=1:size(vardist.means,1)
    for j=1:size(vardist.means,2)
        varmeanNames{i,j} = ['varmean(' num2str(i) ', ' num2str(j) ')'];
    end
    end
    for i=1:size(vardist.means,1)
    for j=1:size(vardist.means,2)
        varcovNames{i,j} = ['varcov(' num2str(i) ', ' num2str(j) ')'];
    end
    end
    names = {varmeanNames{:}, varcovNames{:}}; 
end

% Check if parameters are being optimised in a transformed space.
if ~isempty(vardist.transforms)
  for i = 1:length(vardist.transforms)
    index = vardist.transforms(i).index;
    fhandle = str2func([vardist.transforms(i).type 'Transform']);
    params(index) = fhandle(params(index), 'xtoa');
  end
end
