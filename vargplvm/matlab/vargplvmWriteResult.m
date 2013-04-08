function fileName = vargplvmWriteResult(model, type, dataset, number, varargin)

% SHEFFIELDMLWRITERESULT A simple wrapper for saving an optimised vargplvm model with a consistent filename
% DESC If the first argument is empty, then this function just returns the filename
% SHEFFIELDML 

if ~isfield(model, 'saveName') || isempty(model.saveName)
  if length(dataset) > 0
    dataset(1) = upper(dataset(1));
  end
  if length(type) > 0
    type(1) = upper(type(1));      
  end
  fileName = ['dem' dataset type num2str(number)];
  if nargin > 4
      fileName = [fileName varargin{1}];
  end
  if isoctave
    fileName = [fileName '.mat'];
  end
else
    fileName = model.saveName;
end

  if ~isempty(model)
    save(fileName, 'model');
  end
    
  
