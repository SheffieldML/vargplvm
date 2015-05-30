% VARGPLVMFEATURECLASSIFICATION Use vargplvm features to perform
% discriminative classification with log. regression.
%
%
% ARG: model   : The vargplvm model
% ARG: data    : A struct with:
%                 data.lbls      : Training labels in 1-of-K encoding
%                 data.lblsTs    : Test labels
%                 data.Xpred     : A test latent point (optional if data.Yts given)
%                 data.Xpred_var : A test latent variance (optional if data.Yts given)
%                 data.Yts       : Test output (optional if the 2 above given)
% ARG: options : A struct with options for the function: (Optional)
%                 options.samplesPerObserved: Populate tr. set with this amount of samples per q(x_i)
%                 options.samplesPerOutput  : When predicting FROM q(x*),
%                                             instead predict from many x*_j's sampled from q(x*) and
%                                             average predictions. samplesPerOutput says how many
%                                             x*_j's to sample.
%----
% RET: LogRegError      : Error when not populating training/predictions
% RET: LogRegErrorExt   : Error when populating training set
% RET: LogRegErrorExtOut: Error when populating trainign and predictions
%
% COPYRIGHT: Andreas C. Damianou, 2015
% 
% SEEALSO: vargplvmPredictLatent.m
function [LogRegError, LogRegErrorExt, LogRegErrorExtOut, X_pred, varX_pred] = vargplvmFeatureClassification(model, data, options)

% Train classifiers from X to labels
labels = transformLabels(data.lbls)';
labelsTs = transformLabels(data.lblsTs)';

if ~isfield(data, 'X_pred')
    [X_pred, varX_pred] = ...
        vargplvmPredictLatent(model, data.Yts, [], false, model.globalOpt.reconstrIters,0,[],[],1);
else
    X_pred = data.X_pred;
    varX_pred = data.varX_pred;
end

if nargin < 3 || ~isfield(options, 'samplesPerObserved')
    samplesPerObserved = 5;
else
    samplesPerObserved = options.samplesPerObserved;
end
if nargin < 3 || ~isfield(options, 'samplesPerOutput')
    samplesPerOutput = 10;
else
    samplesPerOutput = options.samplesPerOutput;
end
% Use only SOME dimensions as features (e.g. the ones with best ARD weight)
if nargin < 3 || ~isfield(options, 'dims') || isempty(options.dims)
    dims = 1:model.q;
else
    dims = options.dims;
end

% Populate training set
%--- Take variance into account by sampling new data from the distribution
% Init sizes
Xnew = nan(size(model.vardist.means,1)*samplesPerObserved, size(model.vardist.means,2));
labelsNew = nan(size(Xnew,1), size(labels,2));
k=1;
% Take samples
for n=1:size(model.vardist.means,1)
    for kk = 1:samplesPerObserved
        Xnew(k,:) = model.vardist.means(n,:) + randn(size(model.vardist.means(n,:))).*sqrt(model.vardist.covars(n,:));
        labelsNew(k,:) = labels(n,:);
        k = k + 1;
    end
end
% Augment set with samples
Xext = [model.vardist.means; Xnew];
labelsExt = [labels; labelsNew];
clear 'Xnew' 'labelsNew';

% Training of logistic regression classifier. One for each label
% separately.
nClasses = length(unique(labels));
clear 'B' 'BExt'
for i=1:nClasses
    fprintf('\n # LogReg training for class # %d\n', i)
    lb = zeros(size(model.vardist.means(:,dims),1),1);
    lbExt = zeros(size(Xext(:,dims),1),1);
    lb(labels == i) = 1;
    lbExt(labelsExt == i) = 1;
    B{i} = glmfitWrapper(model.vardist.means(:,dims), lb,'binomial','logit',[],[],[],[],1000);
    BExt{i} = glmfitWrapper(Xext(:,dims), lbExt,'binomial','logit',[],[],[],[],1000);
end

svmmodel    = svmtrain(labels,                 model.vardist.means(:,dims),'-q');
%svmmodelExt = svmtrain(transformLabels(lbExt), model.layer{options.lOut}.vardist.means(:,dims),'-q');

[~, acc,~]    = svmpredict(labelsTs, X_pred(:,dims), svmmodel,'-q');
%[~, accExt] = svmpredict(labelsTs',X_pred, svmmodelExt);

% Prediction of each binary classifier
Ypred_logReg = zeros(size(data.lblsTs));
Ypred_logRegExtOut = zeros(size(data.lblsTs));
Ypred_logRegExt = zeros(size(data.lblsTs));
for i=1:nClasses
    for k=1:samplesPerOutput
        % Sample from the OUTPUT distribution, and then average predictions
        Xsamp = X_pred + randn(size(X_pred)).*sqrt(varX_pred);
        Ypred_logRegExtOut(:,i) = Ypred_logRegExtOut(:,i)+1/samplesPerOutput*glmval(BExt{i}, Xsamp(:,dims), 'logit');
    end
    Ypred_logReg(:,i) = glmval(B{i}, X_pred(:,dims), 'logit')';
    Ypred_logRegExt(:,i) = glmval(BExt{i}, X_pred(:,dims), 'logit')';
end
% Replace predictions with maximum probability (ie, make a decision)
[~,ind]=max(Ypred_logReg');
[~,indExt]=max(Ypred_logRegExt');
[~,indExtOut]=max(Ypred_logRegExtOut');
LogRegError = 0;
LogRegErrorExt = 0;
LogRegErrorExtOut = 0;
for i=1:size(X_pred,1)
    LogRegError = LogRegError + (ind(i) ~= labelsTs(i));
    LogRegErrorExt = LogRegErrorExt + (indExt(i) ~= labelsTs(i));
    LogRegErrorExtOut = LogRegErrorExtOut + (indExtOut(i) ~= labelsTs(i));
end


% fprintf('\n========================== REPORT ===========================\n')
% fprintf('# Error                                             : %d\n', LogRegError);
% fprintf('# Error with %d samp. PerObserved                    : %d\n', samplesPerObserved,LogRegErrorExt);
% fprintf('# Error with %d samp. PerObserved, %d samp.PerOutput : %d\n', samplesPerObserved, samplesPerOutput,LogRegErrorExtOut);
% fprintf('----------------------------------------------------------------\n\n');
% 


N = size(X_pred,1);
fprintf('\n========================== REPORT ===========================\n')
fprintf('# Acc Reg:                                             : %.2f%%\n', (N-LogRegError)/N * 100);
fprintf('# Acc Reg: with %d samp. PerObserved                    : %.2f%%\n', samplesPerObserved,(N-LogRegErrorExt)/N*100);
fprintf('# Acc Reg: with %d samp. PerObserved, %d samp.PerOutput : %.2f%%\n', samplesPerObserved, samplesPerOutput,(N-LogRegErrorExtOut)/N*100);
fprintf('# Acc SVM:                                             : %.2f%%\n', acc(1));
%fprintf('# Acc SVM: with % samp. PerObserved                    : %.2f%%\n', samplesPerObserved,accExt(1));
fprintf('----------------------------------------------------------------\n\n');
