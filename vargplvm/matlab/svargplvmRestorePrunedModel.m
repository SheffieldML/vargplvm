function model = svargplvmRestorePrunedModel(model, Ytr, onlyData, options, parallel)
% SVARGPLVMRESTOREPRUNEDMODEL Restore a pruned shared var-GPLVM model.
% FORMAT
% DESC restores a svargplvm model which has been pruned and it brings it in
% the same state that it was before prunning.
% ARG model: the model to be restored
% ARG Ytr: the training data (it has to be a cell array of length equal to
% model.numModels)
% ARG onlyData: only pruned the data parts. Useful when saving a model which
% is updated after predictions.
% RETURN model : the variational GP-LVM model after being restored
%
% COPYRIGHT: Andreas Damianou,  2011
%
% SEEALSO : svargplvmPruneModel, vargplvmRestorePrunedModel

% VARGPLVM

if nargin <3 || isempty(onlyData)
    onlyData = 0;
end

if nargin <4
    options = [];
end

if nargin < 5 || isempty(parallel)
    parallel = false;
end

try
    pool_open = matlabpool('size')>0;
catch e
    pool_open = 0;
end
if ~pool_open, parallel = false; end

if isfield(model, 'globalOpt')
   if ~iscell(model.globalOpt.balanceModalityDim)
        tmp = repmat(model.globalOpt.balanceModalityDim, 1, model.numModels);
        model.globalOpt.balanceModalityDim = mat2cell(tmp,1,ones(1,size(tmp,2)));
   end
   if isfield(model, 'modalityMapping')
       Ytr = svargplvmMapModalitiesHigh(model.globalOpt.balanceModalityDim, Ytr, model.modalityMapping);
   else
       Ytr = svargplvmMapModalitiesHigh(model.globalOpt.balanceModalityDim, Ytr);
   end
end



if parallel
    modelVardist = model.vardist;
    comp = cell(1,length(model.comp));
    for i=1:model.numModels
        comp{i} = model.comp{i};
    end
    parfor i=1:model.numModels
        comp{i}.vardist = modelVardist;
        comp{i} = vargplvmRestorePrunedModel(comp{i},Ytr{i}, onlyData, options);
    end
    for i=1:model.numModels
        model.comp{i} = comp{i};
    end
    clear 'comp';
else
    if model.numModels > 3 && ~parallel
        pb = myProgressBar(model.numModels, model.numModels);
    end
    for i=1:model.numModels
        model.comp{i}.vardist = model.vardist;
        model.comp{i} = vargplvmRestorePrunedModel(model.comp{i},Ytr{i}, onlyData, options);
        if model.numModels > 3
            pb = myProgressBar(pb,i);
        end
    end
end
