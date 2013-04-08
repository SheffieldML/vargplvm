% Script to load a saved model as well as a predictions file, and
% recalculate the error and restore the prediction matrices.

%% ------- INPUT : settings for *training* the model and load prediction file
% eg.
%{
 dynUsed = 1; % REQUIRED
 experimentNo = 10001;
 dataSetName = 'missa';
 dataSetSplit = 'blocks';
 blockSize = 5;
 cutP = [35:165 250:288];
 % Load prediction file
 load demMissaVargplvm1013Pred;
%}
%{
dynUsed=1;
experimentNo = 10018;
dataSetName = 'missa';
dataSetSplit = 'blocks';
blockSize = 5; 
cutP = [35:165 250:288];
load demMissaVargplvm10018Pred;
%}
%% -----  SCRIPT -----------------------------------------------
trainModel=false;
predWithMs = false;
fileToSave = [];

if dynUsed
    demHighDimVargplvm3
else
    demHighDimVargplvmStatic
end


%--- Run the part which loads and prepares the test data

    display = 1;
    indexP = [];
    Init = [];
    
    fprintf(1, '# Partial reconstruction of test points...\n');
    
    demHighDimPrepareTestData
    
    Yts(:,indexMissing) = NaN;


%---- Restore modelUpdated and do predictions
if dynUsed
    modelUpdated = vargplvmRestorePrunedModel(prunedModelUpdated, Ytr,true);
    Varmu = vargplvmPosteriorMeanVar(modelUpdated, Testmeans, Testcovars);
else
    Varmu = vargplvmPosteriorMeanVar(model, Testmeans, Testcovars);
end


% Mean error per pixel
errorFull = sum(sum( abs(Varmu(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
fprintf(1,'# GPLVM Error (in the missing dims) with missing inputs:%d\n', errorFull);

if dynUsed
    errOnlyTimes2 = sum(sum( abs(Varmu2(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
    fprintf(1,'# GPLVM Error (in the missing dims) only times:%d\n', errOnlyTimes2);
end

errorFullPr = sum(sum( abs(Varmu(:,indexPresent) - YtsOriginal(:,indexPresent)) ))/prod(size(YtsOriginal(:,indexPresent)));
fprintf(1,'# GPLVM Error (in the present dims) with missing inputs:%d\n', errorFullPr);

if dynUsed
    errorFullPr2 = sum(sum( abs(Varmu2(:,indexPresent) - YtsOriginal(:,indexPresent)) ))/prod(size(YtsOriginal(:,indexPresent)));
    fprintf(1,'# GPLVM Error (in the present dims) only times:%d\n', errorFullPr2);
end


if exist('skipNN') && skipNN
    return
end

%--- -------- NN  ----------
meanFrame = repmat(mean(Ytr),size(YtsOriginal,1),1);
errMeanFrame = sum(sum(abs(meanFrame(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
fprintf(1,'# Mean Frame Error (in the missing dims) with missing inputs: %f\n', errMeanFrame);

%{
fprintf('# Sequential NN...\n')
try
    [NNpredSeq, errorsNNSeq] = NNseq(Ytr, YtsOriginal, Yts,1:6);
    [minEr, indMinEr] = min(errorsNNSeq);
    fprintf(1,'# NN Seq.(%d) Error (in the missing dims) with missing inputs:%f\n', indMinEr, minEr);
catch e
end
%}
fprintf(1,'# Classic NN prediction...\n');

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
if ~timeNN
    bestK=2; bestErr = errorNNPart; NNmuPartBest = NNmuPart;
    clear NNmuPart %%%
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

%% Movies
% playMov(h,w,0.5,Varmu, NNmuPartBest )