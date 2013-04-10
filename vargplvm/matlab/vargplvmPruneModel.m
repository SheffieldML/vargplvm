function mm = vargplvmPruneModel(model, onlyData)
% VARGPLVMPRUNEMODEL Prune a var-GPLVM model.
% FORMAT
% DESC prunes a VAR-GPLVM model by removing some fields which can later be
% reconstructed based on what is being kept. Used when storing a model.
% ARG model: the model to be pruned
% ARG onlyData: only prune the data parts. Useful when saving a model which
% is updated after predictions.
% RETURN mm : the variational GP-LVM model after being pruned
%
% COPYRIGHT: Andreas Damianou, Michalis Titsias, Neil Lawrence, 2011
%
% SEEALSO : vargplvmReduceModel, vargplvmRestorePrunedModel

% VARGPLVM

if exist('onlyData') && onlyData
    model.m = [];
    model.y = [];
   if isfield(model, 'mOrig')
        model.mOrig = [];
   end
   if isfield(model, 'scale')
       model.scale = [];
   end
   if isfield(model, 'bias')
       model.bias = [];
   end
   model.P = [];
   model.B = [];
   mm = model;
   return
end


fieldsToKeep = ...
    {'type','approx','learnScales', 'optimiseBeta','betaTransform','q','d','N','optimiser','nPrivateParams','learnSigmaf',...
    'bias','scale','DgtN', 'kern','k','fixInducing','inducingIndices','X_u','beta','prior','vardist','numParams', ...
    'nParams','learnBeta', 'dynamics', 'iters', 'date', 'fixedBetaIters', 'dataSetInfo','id', ...
    'initVardist', 'initVardistIters', 'saveName', 'globalOpt', ...
    'latentIndices'}'; % This line for hsvargplvm

mm=[];
for i=1:length(fieldsToKeep)
    if isfield(model, fieldsToKeep{i})
        f = getfield(model,fieldsToKeep{i});
        mm=setfield(mm,fieldsToKeep{i},f);
    else
        %if ~strcmp(fieldsToKeep{i}, 'inducingIndices') && ~strcmp(fieldsToKeep{i}, 'iters') && ~strcmp(fieldsToKeep{i}, 'fixedBetaIters')
        %    fprintf(['??? Field ' fieldsToKeep{i} ' was missing from model! This field was skipped...\n']);
        %end
    end
end

if isfield(model, 'dynamics') && ~isempty(model.dynamics)
    mm.dynamics.X = [];
end