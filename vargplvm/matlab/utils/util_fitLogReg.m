function [LogRegError, Ypred_logReg, B, ind] = util_fitLogReg(curInp, curOut, curInpTs, lblsTs, pcaDim, alwaysPCA)

if nargin < 3, curInpTs = []; end
if isempty(curInpTs)
    assert(nargin < 4 || isempty(lblsTs));
end
if nargin < 5 || isempty(pcaDim), pcaDim = size(curInp,1)-1; end
if nargin < 6 || isempty(alwaysPCA), alwaysPCA = false; end

nClasses = length(unique(curOut));
if size(curInp,1)<size(curInp,2) || alwaysPCA
    warning('Low rank matrix! Applying PCA on the inputs...')
    [~,V] = pca(curInp, pcaDim);
    curInp = curInp * V;
    if ~isempty(curInpTs)
        curInpTs = curInpTs * V;
    end
end
for i=1:nClasses
    fprintf('\n # LogReg training for class # %d\n', i)
    lb = zeros(size(curInp,1),1);
    lb(curOut == i) = 1;
    B{i} = glmfit(curInp, lb,'binomial','logit'); % Logistic regression
end
if ~isempty(curInpTs)
    Nts = size(lblsTs,1);
    labelsTs = transformLabels(lblsTs);
    % Prediction of each binary classifier
    Ypred_logReg = zeros(size(lblsTs));
    for i=1:nClasses
        Ypred_logReg(:,i) = glmval(B{i},curInpTs,'logit')';
    end
    % Replace predictions with maximum probability (ie, make a decision)
    [~,ind]=max(Ypred_logReg');
    LogRegError = 0;
    for i=1:Nts
        LogRegError = LogRegError + (ind(i) ~= labelsTs(i));
    end
else
    LogRegError = [];
end