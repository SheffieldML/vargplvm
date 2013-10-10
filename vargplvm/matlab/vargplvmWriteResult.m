function fileName = vargplvmWriteResult(model, type, dataset, number, varargin)

% VARGPLVMWRITERESULT A simple wrapper for saving an optimised vargplvm model with a consistent filename
% DESC If the first argument is empty, then this function just returns the filename
% VARGPLVM 

if ~isfield(model, 'saveName') || isempty(model.saveName)
  if length(dataset) > 0
    dataset(1) = upper(dataset(1));
  end
  if length(type) > 0
    type(1) = upper(type(1));      
  end
  fileName = ['dem' dataset type num2str(number)];
  if nargin > 4
      fileName = [fileName varargin{1}]; % Potential suffix
  end
  if nargin > 5
      fileName = [varargin{2} fileName]; % Potential prefix
  end
  if isoctave
    fileName = [fileName '.mat'];
  end
else
    fileName = model.saveName;
end

if ~isempty(model)
    try
        save(fileName, 'model');
        fprintf('# Saved ''%s.mat'' !\n', fileName)
    catch e
        e.getReport
        if ~exist('./tmp.mat', 'file')
            save('tmp.mat', 'model');
            warning(0,['Could not save ' fileName '! Saved ''tmp.mat'' instead!'])
        end
    end
end


