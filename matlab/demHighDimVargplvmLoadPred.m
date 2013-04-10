% DEMHIGHDIMVARGPLVMLOADPRED Load predictions of the var. GPLVM made on high dimensional video datasets.
% DESC: As the dataset is high dimensional it is inefficient to store the predictions made by the model.
% Instead, we store minimal information and recompute the predictions very fast, without needing to do
% the optimisation again. The xxxxPred file which must be loaded must contain the following quantities:
% Testmeans, Testcovars, prunedModelUpdated (see demHighDimVargplvm3.m for details). These files are
% related to q(X_*).
% This script then does the necessary calculations (restores the model and calculates the posterior)
% to find the predictions.
%
% COPYRIGHT: Andreas C. Damianou, Michalis K. Titsias, 2010 - 2011
% SEEALSO: demHighDimVargplvm3.m
% VARGPLVM

%---- These two steps must be done a priori ----
% <dataset info and experimentNo>:
% 	Here, set the dataSetName, experimentNo, and any other options
% 	related to the dataset.
% load xxxxPred % load the .mat file with the predictions
%-------------

predWithMs=0; doSampling=0;
demHighDimVargplvmTrained

model.y = Ytr;
model.m= gpComputeM(model); %%%
model.y=[];



Nstar = size(YtsOriginal,1);
modelUpdated = vargplvmRestorePrunedModel(prunedModelUpdated,Ytr,1);
modelUpdated.m = modelUpdated.mOrig;
modelUpdated.P = modelUpdated.P1 * (modelUpdated.Psi1' * modelUpdated.m);
[Testmeans2 Testcovars2] = vargplvmPredictPoint(model.dynamics, timeStampsTest);
Varmu2 = vargplvmPosteriorMeanVar(model, Testmeans2, Testcovars2);


display = 1;
indexP = [];
Init = [];

% w=360; % width
% h=288; % height
w=width; h=height;


if ~exist('cut')
    cut = 'vertically';
    if strcmp(dataSetName,'missa') || strcmp(dataSetName,'ADN') || strcmp(dataSetName,'claire')
        cut = 'horizontally';
    end
end

switch cut
    case 'horizontally'
        if strcmp(dataSetName,'ADN')
            cutPoint=round(h/1.55);
        elseif strcmp(dataSetName,'claire')
            cutPoint=round(h/2);
        elseif strcmp(dataSetName, 'grandma')
            cutPoint=round(h/1.68);
        else %missa
            cutPoint=round(h/2)+13;
        end
        if exist('cutP')    cutPoint = cutP; end
        mask = [ones(1,cutPoint) zeros(1,h-cutPoint)];
        mask=repmat(mask, 1,w);
        indexMissing = find(mask);
    case 'vertically'
        if strcmp(dataSetName,'missa')
            indexMissing=1:round(size(Yts,2)/1.8);
        elseif strcmp(dataSetName,'dog')
            indexMissing = 1:round(size(Yts,2)/1.70);
        elseif strcmp(dataSetName,'ADN')
            indexMissing = 1:round(size(Yts,2)/1.7);
        elseif strcmp(dataSetName,'ocean')
            indexMissing = 1:round(size(Yts,2)/1.6);
        elseif strcmp(dataSetName,'head')
            indexMissing=1:round(size(Yts,2)/2.08);
        else
            indexMissing=1:round(size(Yts,2)/2);
        end
end
indexPresent = setdiff(1:model.d, indexMissing);

%
Yts = YtsOriginal;
Yts(:,indexMissing) = NaN;

Varmu = vargplvmPosteriorMeanVar(modelUpdated, Testmeans, Testcovars);
% Mean error per pixel
errorFull = sum(sum( abs(Varmu(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
fprintf(1,'# GPLVM Error (in the missing dims) with missing inputs:%d\n', errorFull);

errOnlyTimes2 = sum(sum( abs(Varmu2(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
fprintf(1,'# GPLVM Error (in the missing dims) only times:%d\n', errOnlyTimes2);

errorFullPr = sum(sum( abs(Varmu(:,indexPresent) - YtsOriginal(:,indexPresent)) ))/prod(size(YtsOriginal(:,indexPresent)));
fprintf(1,'# GPLVM Error (in the present dims) with missing inputs:%d\n', errorFullPr);

errorFullPr2 = sum(sum( abs(Varmu2(:,indexPresent) - YtsOriginal(:,indexPresent)) ))/prod(size(YtsOriginal(:,indexPresent)));
fprintf(1,'# GPLVM Error (in the present dims) only times:%d\n', errorFullPr2);


%-------- NN  ----------
fprintf(1,'# NN(2) prediction...\n');
mini =[];
sortedInd = zeros(size(Yts,1), size(Ytr,1));
for i=1:size(Yts,1)
    dst = dist2(Yts(i,indexPresent), Ytr(:,indexPresent));
    % mind: smaller distance %mini(i): smaller index
    %[mind, mini(i)] = min(dst);
    [void ,sortedInd(i,:)] = sort(dst,2);
end
clear dst %%%


NNmuPart = zeros(Nstar, size(Ytr,2));
% Set to 1 to find the two NN's in time space. set to 0 to find in
% data space. timeNN=1 should be set only for datasets created with
% the "everyTwo" option for the split.
timeNN=0;
k=2; % the k parameter in k-NN
for i=1:Nstar
    if timeNN
        if i < Nstar
            NNmuPart(i, indexMissing) = 0.5*(Ytr(i,indexMissing) + Ytr(i+1,indexMissing));
        else
            NNmuPart(i, indexMissing) = Ytr(end,indexMissing);
        end
    else
        NNmuPart(i,indexMissing) = Ytr(sortedInd(i,1),indexMissing);
        for n=2:k
            NNmuPart(i,indexMissing) = NNmuPart(i,indexMissing)+Ytr(sortedInd(i,n),indexMissing);
        end
        NNmuPart(i,indexMissing) = NNmuPart(i,indexMissing)*(1/k);
    end
    NNmuPart(i, indexPresent) = Yts(i, indexPresent);
end

% Mean absolute error per pixel
%errorNNPart = mean(abs(NNmuPart(:) - YtsOriginal(:)));
errorNNPart = sum(sum( abs(NNmuPart(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
fprintf(1,'# NN(2) Error (in the missing dims) with missing inputs:%d\n', errorNNPart);


%%%%%%  Try different k for k-NN to see which is best
fprintf(1,'# k-NN prediction for k=[1 3 4 5]...\n');
if ~timeNN
    bestK=2; bestErr = errorNNPart; NNmuPartBest = NNmuPart;
    for k=[1 3 4 5]; % the k parameter in k-NN, try different values
        NNmuPartK = zeros(Nstar, size(Ytr,2));
        for i=1:Nstar
            NNmuPartK(i,indexMissing) = Ytr(sortedInd(i,1),indexMissing); % first NN
            for n=2:k % more NN's, if k>1
                NNmuPartK(i,indexMissing) = NNmuPartK(i,indexMissing)+Ytr(sortedInd(i,n),indexMissing);
            end
            NNmuPartK(i,indexMissing) = NNmuPartK(i,indexMissing)*(1/k); % normalise with the number of NN's
            
            NNmuPartK(i, indexPresent) = Yts(i, indexPresent);
        end
        errorNNPartK = sum(sum( abs(NNmuPartK(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
        if errorNNPartK < bestErr
            bestErr = errorNNPartK;
            bestK=k;
            NNmuPartBest = NNmuPartK;
        end
    end
    clear NNmuPartK
end
% Mean absolute error per pixel
fprintf(1,'# NNbest(%d) Error (in the missing dims) with missing inputs:%d\n',bestK, bestErr);
%%%%%%%%%%%%%
