function model = svargplvmRestorePrunedModel(model, Ytr, onlyData, options)
% SHEFFIELDMLRESTOREPRUNEDMODEL Restore a pruned shared var-GPLVM model.
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

% SHEFFIELDML

if nargin <3
    onlyData = 0;
end

if nargin <4
    options = [];
end

for i=1:model.numModels
    model.comp{i}.vardist = model.vardist;
    model.comp{i} = vargplvmRestorePrunedModel(model.comp{i},Ytr{i}, onlyData, options);
end

