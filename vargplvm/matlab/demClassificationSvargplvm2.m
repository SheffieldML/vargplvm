% DEMCLASSIFICATIONSVARGPLVM2
% COPYRIGHT Andreas C. Damianou, 2012
% SEEALSO: demClassificationSvargplvm, demClassificationSvargplvm3,  demClassificationSvargplvm4
% VARGPLVM

% Like demClassificationSvargplvm but we train
% with RANDOM subsets of the training data and in the end we take the mean.

% Fix seeds

if ~exist('onlyTest'), onlyTest = false; end
if ~exist('doPredictions'),    doPredictions = false; end
if ~exist('testDataPerClass'), testDataPerClass = -1; end %%% All data
if ~exist('printPlots'), printPlots = true; end

dataPerClass = -1; %%%% ALL data

addpath(genpath('../../bc-vargplvm/matlab'));
addpath(genpath('../../vargplvm/matlab'));

% This script initialises the options structure 'globalOpt'.
svargplvm_init;

if onlyTest
    load(['demOilSvargplvm' num2str(experimentNo)]);
    globalOpt = model.globalOpt;
end


globalOpt.dataSetNames = {'oilData', 'oilLabels'};
globalOpt.dataSetName = 'oil';
globalOpt.dataToKeep = -1; %%


if exist('timeStamps')
    [globalOpt, YtrAll] = bc_loadData(globalOpt, timeStamps);
else
    [globalOpt, YtrAll] = bc_loadData(globalOpt);
end
Yall{1} = YtrAll.Y;
lblsTrain = YtrAll.lbls;
labelsTrain = transformLabels(lblsTrain);
% ---- We want labels with -1 and 1, not 0 and 1 (anti-correlation)
YtrAll.lbls(YtrAll.lbls == 0) = -1;
Yall{2} = YtrAll.lbls;



globalOptTest = globalOpt;
globalOptTest.dataPerClass = testDataPerClass;
globalOptTest.dataSetName = 'oilTest';
if exist('timeStamps')
    [globalOptTest, YtsAll] = bc_loadData(globalOptTest, timeStamps);
else
    [globalOptTest, YtsAll] = bc_loadData(globalOptTest);
end
lblsTest = YtsAll.lbls;
labelsTest = transformLabels(lblsTest);
Yts{1} = YtsAll.Y;
YtsAll.lbls(YtsAll.lbls == 0) = -1;
Yts{2} = YtsAll.lbls;


if isfield(globalOpt, 'normaliseData') && globalOpt.normaliseData
    Yall{1} = utils_normaliseData(Yall{1});
    Yts{1} = utils_normaliseData(Yts{1});
end


%%%%%%%%%%
indTemp = randperm(size(Yall{1},1));
indTrial = indTemp(1:dataToKeep);
Yall{1} = Yall{1}(indTrial,:);
Yall{2} = Yall{2}(indTrial,:);
lblsTrain = lblsTrain(indTrial,:);
labelsTrain = labelsTrain(indTrial);
%%%%%%%%%%%



if onlyTest
    model = svargplvmRestorePrunedModel(model, Yall);
end

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
        for i=1:numberOfDatasets
            % Set up dynamcis (i.e. back-constraints) model
            optionsDyn{i}.type = 'vargpTime';
            optionsDyn{i}.inverseWidth=30;
            %   optionsDyn.vardistCovars = vardistCovarsMult;
            optionsDyn{i}.initX = globalOpt.initX;
            optionsDyn{i}.constrainType = globalOpt.dynamicsConstrainType;
            
            if exist('timeStampsTraining')
                optionsDyn{i}.t = timeStampsTraining;
            end
            if exist('labelsTrain') && ~isempty(labelsTrain)
                optionsDyn{i}.labels = labelsTrain;
            end
        end
    else
        optionsDyn= [];
    end
    
    
    
    
    model = svargplvmModelCreate(Ytr, globalOpt, options, optionsDyn);
    if exist('diaryFile')
        model.diaryFile = diaryFile;
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
        fprintf('aaa\n\n\n')
        model.dynamics.kern.comp{2}.variance = whiteVar;
    end
    %%%%
    
    
    % Force kernel computations
    params = svargplvmExtractParam(model);
    model = svargplvmExpandParam(model, params);
    

    
    %%
    fprintf('# Median of vardist. covars: %d \n',median(median(model.vardist.covars)));
    fprintf('# Min of vardist. covars: %d \n',min(min(model.vardist.covars)));
    fprintf('# Max of vardist. covars: %d \n',max(max(model.vardist.covars)));
    
    
    
    model = svargplvmOptimiseModel(model);
    
    svargplvmShowScales(model)
    figure
end

%%
    capName = model.globalOpt.dataSetName;
    capName(1) = upper(capName(1));
    errors = fgplvmNearestNeighbour(model.comp{1}, lblsTrain);
    model2 = vargplvmReduceModel(model.comp{1}, 2);
    errors2 = fgplvmNearestNeighbour(model2, lblsTrain);
    fprintf('# Visualisation errors in all dims/in 2D:  %d / %d\n', errors, errors2)
    if printPlots    
        vargplvmPrintPlot(model2, lblsTrain, [capName 'Vargplvm'], model.globalOpt.experimentNo);
               
        %%
        % For 3D plots:
        % labels = model.options{1}.labels; dims = [1 3 5];
        % plot3k({model.X(:,dims(1)) model.X(:,dims(2)) model.X(:,dims(3))}, 'ColorData', labels, 'Marker', {'x',6});
    end

%%
obsMod = 1; % one of the involved sub-models (the one for which we have the data)
infMod = setdiff(1:2, obsMod);

sharedDims = svargplvmFindSharedDims(model);

%---------------------------- PREDICTIONS ---------------
if ~doPredictions
    return
end


%%

if ~exist('testOnTraining')
    testOnTraining=0;
end


%--------------------%
svargplvmPredictions %---- Script returning: ZpredMuAll and mini(the indices for NN)
%--------------------%


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
