% DEMCLASSIFICATIONGENERALSVARGPLVM
% COPYRIGHT Andreas C. Damianou, 2012
% SEEALSO: demClassificationSvargplvm, demClassificationSvargplvm2, demClassificationSvargplvm3,  demClassificationSvargplvm4
% VARGPLVM

%-- Example script
%{
clear
dynUsed = 0;
dataToKeep = 60;
dataSetName = 'mulan_emotions';
baseKern = {{'linard2', 'white','bias'},{'linard2', 'white','bias'}}; % {'rbfardjit', 'rbfardjit'}; 
indPoints = 120;
initX = 'ppca';
initial_X = 'concatenated';
% initial_X = 'separately'; latentDimPerModel = 6;
initSNR = 100;
latentDim = 15;
initVardistIters = 60;
itNo = 80;
scale2var1 = 1;
doPredictions = true;
dynamicsConstrainType = {}; %%% dynUsed = 0

%---
if dynUsed
    dynamicsConstrainType = {'labels'}; % dynIsed = 1
    vardistCovarsMult = 0.7;
    dynamicKern = {'lin', 'white', 'bias'};
end
%}
%------

%% 

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);
if ~exist('onlyTest'), onlyTest = false; end
if exist('diaryFile'),    diary(diaryFile), end
if ~exist('doPredictions'),    doPredictions = false; end
if ~exist('testDataPerClass'), testDataPerClass = -1; end
if ~exist('displayOutput'), displayOutput = true; end
if ~exist('testDataToKeep'), testDataToKeep = -1; end


addpath(genpath('../../vargplvm/matlab'));

% This script initialises the options structure 'globalOpt'.
svargplvm_init;

if onlyTest
    fName = vargplvmWriteResult([], 'svargplvm',dataSetName,experimentNo);
    load(fName);
    globalOpt = model.globalOpt;
    clear model
end

globalOpt.dataSetNames = {[globalOpt.dataSetName 'Data'], [globalOpt.dataSetName 'Labels']};


[Yall{1}, Yall{2}, Yts{1}, Yts{2}] = svargplvmLoadData(globalOpt.dataSetName);
%TEMP
if exist('splitTraining') && splitTraining
    th = 200;
    Yts{1} = Yall{1}(th+1:end,:);
    Yts{2} = Yall{2}(th+1:end,:);
    Yall{1} = Yall{1}(1:th,:);
    Yall{2} = Yall{2}(1:th,:);
end

if testDataToKeep ~= -1
    inds = randperm(size(Yts{2},1));
    Yts{1} = Yts{1}(inds(1:testDataToKeep),:);
    Yts{2} = Yts{2}(inds(1:testDataToKeep),:);
end


%---
if globalOpt.dataToKeep ~= -1
    Yall{1} = Yall{1}(1:globalOpt.dataToKeep, :);
    Yall{2} = Yall{2}(1:globalOpt.dataToKeep, :);
end

if isfield(globalOpt, 'transformMultiLabel') && globalOpt.transformMultiLabel
    [Yall{1}, Yall{2}] = utils_transformMultiLabel(Yall{1}, Yall{2});
end

%---
%-- Sort the data by transforming each row to number)
[void, indexSort] = sort(bin2dec(num2str(Yall{2})));
Yall{1} = Yall{1}(indexSort,:);
Yall{2} = Yall{2}(indexSort,:);
%Yall{2} = Yall{1};%%%% TEMP %%%%%%%%%%%%%%
[void, indexSort] = sort(bin2dec(num2str(Yts{2})));
Yts{1} = Yts{1}(indexSort,:);
Yts{2} = Yts{2}(indexSort,:);


%--- if we want labels with -1 and 1 instead of 0 and 1 (forces anti-correlation)
if exist('kernelAnticorrelation') && kernelAnticorrelation
    Yall{2}(Yall{2} == 0) = -1;
    Yts{2}(Yts{2} == 0) = -1;
end
lblsTrain = Yall{2};
lblsTest = Yts{2};
labelsTrain = transformLabels(lblsTrain);
labelsTest = transformLabels(lblsTest);



numberOfDatasets = length(Yall);

% globalOpt.baseKern = {'rbfardjit', 'linard2'};

%globalOpt.baseKern = {'rbfardjit', 'rbfardjit'};

globalOpt.indPoints = min(globalOpt.indPoints, size(Yall{1},1));


%-- Load datasets
for i=1:numberOfDatasets
    Y = Yall{i};
    dims{i} = size(Y,2);
    N{i} = size(Y,1);
    indTr = globalOpt.indTr;
    if indTr == -1
        indTr = 1:N{i};
    end
    if ~exist('Yts')
        indTs = setdiff(1:size(Y,1), indTr);
        Yts{i} = Y(indTs,:);
    end
    Ytr{i} = Y(indTr,:);
    
    t{i} = linspace(0, 2*pi, size(Y, 1)+1)'; t{i} = t{i}(1:end-1, 1);
    timeStampsTraining{i} = t{i}(indTr,1); %timeStampsTest = t(indTs,1);
end

for i=2:numberOfDatasets
    if N{i} ~= N{i-1}
        error('The number of observations in each dataset must be the same!');
    end
end

%%
if ~onlyTest
    % Free up some memory
    clear('Y')
    
    options = svargplvmOptions(Ytr, globalOpt, labelsTrain);
    
    
    
    if ~isempty(globalOpt.dynamicsConstrainType)
        dynUsed=true;
        % Set up dynamcis (i.e. back-constraints) model
        optionsDyn.type = 'vargpTime';
        optionsDyn.inverseWidth=30;
        %   optionsDyn.vardistCovars = vardistCovarsMult;
        optionsDyn.initX = globalOpt.initX;
        optionsDyn.constrainType = globalOpt.dynamicsConstrainType;
        
        if exist('timeStampsTraining')
            optionsDyn.t = timeStampsTraining;
        end
        if exist('labelsTrain') && ~isempty(labelsTrain)
            optionsDyn.labels = labelsTrain;
        end
    else
        dynUsed = false;
        optionsDyn= [];
    end
    
    model = svargplvmModelCreate(Ytr, globalOpt, options, optionsDyn);
    if exist('diaryFile')
        model.diaryFile = diaryFile;
    end
    
    if globalOpt.equalizeScales %% if not needed
        model = svargplvmEqualizeScales(model, globalOpt);
    end 

    
    model.globalOpt = globalOpt;
    model.options = options;
    
    
    
    %-- Define what level of parallelism to use (w.r.t submodels or/and w.r.t
    % datapoints).
    %{
fprintf('# Parallel computations w.r.t the submodels!\n');
model.parallel = 1;
model = svargplvmPropagateField(model,'parallel', 1);
%
fprintf('# Parallel computations w.r.t the datapoints!\n');
model.vardist.parallel = 1;
for i=1:model.numModels
    model.comp{i}.vardist.parallel = 1;
end
    %}
    
            %%%% TEMP
    if exist('whiteVar')
        model.dynamics.kern.comp{2}.variance = whiteVar;
    end
    %%%%
    
    
    % Force kernel computations
    params = svargplvmExtractParam(model);
    model = svargplvmExpandParam(model, params);
    

    
    %%
    if dynUsed
        fprintf('# Median of vardist. covars: %d \n',median(median(model.vardist.covars)));
        fprintf('# Min of vardist. covars: %d \n',min(min(model.vardist.covars)));
        fprintf('# Max of vardist. covars: %d \n',max(max(model.vardist.covars)));
    end
    
    
    if displayOutput
        model = svargplvmOptimiseModel(model);
        figure
    else
        model = svargplvmOptimiseModelNoDisplay(model);
    end
    svargplvmShowScales(model, displayOutput)
    sharedDims = svargplvmFindSharedDims(model)
    
else
    load(fName);
    model = svargplvmRestorePrunedModel(model, Yall);
end



%{
%     capName = model.globalOpt.dataSetName;
%     capName(1) = upper(capName(1));
%     errors = fgplvmNearestNeighbour(model.comp{1}, lblsTrain);
%     model2 = vargplvmReduceModel(model.comp{1}, 2);
%     errors2 = fgplvmNearestNeighbour(model2, lblsTrain);
%     fprintf('# Visualisation errors in all dims/in 2D:  %d / %d\n', errors, errors2)
%     if displayOutput    
%         vargplvmPrintPlot(model2, lblsTrain, [capName 'Vargplvm'], model.globalOpt.experimentNo);
%                
%         %%
%         % For 3D plots:
%          labels = labelsTrain; dims = [1 2 3];
%          plot3k({model.X(:,dims(1)) model.X(:,dims(2)) model.X(:,dims(3))}, 'ColorData', labels, 'Marker', {'x',6});
%     end
%}
%%
obsMod = 1; % one of the involved sub-models (the one for which we have the data)
infMod = setdiff(1:2, obsMod);

sharedDims = svargplvmFindSharedDims(model);

%---------------------------- PREDICTIONS ---------------
if ~doPredictions
    return
end


%%

%---
if isfield(model, 'dynamics')
    modelOld = model;
    model.dynamicsRemoved = model.dynamics; %%%%%%%%5
    model = rmfield(model, 'dynamics');
    for j=1:model.numModels
        model.comp{j} = rmfield(model.comp{j}, 'dynamics');
    end
    
    if exist('reoptimiseInTest') && ~isempty(reoptimiseInTest)
        model = svargplvmReOptimiseModel(model, reoptimiseInTest{1}, reoptimiseInTest{2}, experimentNo*10);
		svargplvmShowScales(model,false)
		sharedDims = svargplvmFindSharedDims(model);
    end
    
end
%---

%--------------------%
svargplvmPredictions2 %---- Script returning: ZpredMuAll and mini(the indices for NN)
%--------------------%

%---
if exist('splitTraining') && splitTraining
    [bestThreshold, bestError, totalError] = svargplvmOptimiseThreshold(ZpredMuAllOrig, Yts{2});
    save(['optThreshold' num2str(experimentNo/1000) '.mat'], 'bestThreshold', 'bestError', 'totalError');
    fprintf('Optimised threshold error: %.5f\n', totalError);
else
    try
        load(['optThreshold' num2str(experimentNo)])
        ss = 0;
        ZpredAll = [];
        for d=1:size(ZpredMuAllOrig,2)
            Ztrue = Yts{2}(:,d);
            Zpred = ZpredMuAllOrig(:,d);
            Zpred(Zpred >= bestThreshold(d)) = 1;
            Zpred(Zpred < bestThreshold(d)) = 0;
            ZpredAll = [ZpredAll Zpred];
            hammingLossOpt = sum(sum(abs(Zpred - Ztrue))) / (size(Zpred,1) * size(Ztrue,2));
            ss = ss+hammingLossOpt;
        end
        hammingLossOpt = ss / size(Yts{2},2);
        fprintf('Error with optimal threshold: %.5f\n', hammingLossOpt)
    catch e
    end
end
%---


    hammingLoss = sum(sum(abs(ZpredMuAll - Yts{2}))) / (size(ZpredMuAll,1) * size(Yts{2},2));
    fprintf('Error: %.5f\n', hammingLoss)
    
    try
        predictions.ZOrig = ZpredMuAllOrig;
        predictions.Z = ZpredMuAll;
        predictions.X = XpredAll;
        predictions.varZ = varZpredAll;
        save(['demMulan_emotionsSvargplvmPred' num2str(experimentNo) '.mat'], 'hammingLoss', 'predictions');
    catch e
        getReport(e)
    end

%{
%--- For classification we are interested in labels
ZpredMuAll(ZpredMuAll > 0) = 1;
ZpredMuAll(ZpredMuAll < 0) = -1;
NNpred = Ytr{infMod}(mini,:);
realLabels = lblsTest(testInd,:);

% Find the error in the labels ------------
[labelErrors, labelSharingErrors, wrongIndices] = findLabelErrors(realLabels, ZpredMuAll);
[labelErrorsNN, void, wrongIndicesNN] = findLabelErrors(realLabels,NNpred);

fprintf('# Gplvm label errors: %d\n', labelErrors);
fprintf('# Gplvm label sharing errors: %d\n', labelSharingErrors);
fprintf('# NN label errors:    %d\n', labelErrorsNN);

% Confusion matrix
cmat_a = confMatrix( XpredAll, model.X, labelsTrain, labelsTest, 1);
cmat_rel = cmat_a ./ repmat(sum(cmat_a, 2), 1, length(unique(labelsTrain)));
cmat_rel
%}
