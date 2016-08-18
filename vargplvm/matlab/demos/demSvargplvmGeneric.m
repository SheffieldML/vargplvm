% DEMSVARGPLVMGENERIC A generic demo for the MRD method.
%
% COPYRIGHT: Andreas C. Damianou, 2012
%
% SEEALO: MRDEmbed.m
%
% VARGPLVM

%------
% This is a generic demo for svargplvm. It only requires a cell-array
% Yall to be a priori created. Each cell should contain one dataset.
% This demo will use default options and configurations for the model.
% These can be changed by investigating the fields of globalOpt structure,
% which are set during 'svargplvm_init'.

%{ 
%eg:
clear
dataSetNames = 'toy'; toyDataCreate = 'fols'; svargplvmPrepareData
baseKern = {{'linard2','white'},{'linard2','white'}};
latentDimPerModel = 4;
init_X = 'separately';
initVardistIters = 200;
itNo = 200;
demSvargplvmGeneric
%}

% % Fix seeds
% randn('seed', 1e5);
% rand('seed', 1e5);

if ~exist('Yall', 'var')
    error('# This demo requires that a cell-array Yall is a priori created. This array should contain one dataset per cell.')
end
if ~exist('trainModel', 'var'), trainModel = true; end

M = length(Yall);

% This script initialises the options structure 'globalOpt'.
svargplvm_init;

% Don't use more inducing points than data
globalOpt.indPoints = min(globalOpt.indPoints, size(Yall{1},1));

%-- Load datasets
for i=1:M
    Y = Yall{i};
    dims{i} = size(Y,2);
    if i == 1
        Nall = size(Y,1);
    else
        if size(Y,1) ~= Nall
            error('The number of observations in each dataset must be the same!');
        end
    end
    indTr = globalOpt.indTr;
    if indTr == -1
        indTr = 1:Nall;
    end
    if ~exist('Yts')
        indTs = setdiff(1:size(Y,1), indTr);
        Yts{i} = Y(indTs,:);
    end
    Ytr{i} = Y(indTr,:);
end
clear('Y','Yall')

%--- Create model
if ~exist('options','var')
    options = svargplvmOptions(Ytr, globalOpt);
else
    warning('options field exists already!');
end

if ~isempty(globalOpt.dynamicsConstrainType)
    if exist('optionsDyn', 'var')
        warning('optionsDyn field exists already!');
    else
        if exist('timeStampsTraining', 'var')
            t = timeStampsTraining;
        else
            t = linspace(0, 2*pi, Nall+1)'; t = t(1:end-1, 1);
        end
        timeStampsTraining = t(indTr,1); %timeStampsTest = t(indTs,1);
        optionsDyn.type = 'vargpTime';
        optionsDyn.inverseWidth=30;
        optionsDyn.initX = globalOpt.initX;
        optionsDyn.constrainType = globalOpt.dynamicsConstrainType;
        if exist('timeStampsTraining', 'var')
            optionsDyn.t = timeStampsTraining;
        end
        if exist('labelsTrain', 'var') && ~isempty(labelsTrain)
            optionsDyn.labels = labelsTrain;
        end
    end
else
    optionsDyn= [];
end

model = svargplvmModelCreate(Ytr, globalOpt, options, optionsDyn);
if exist('diaryFile'), model.diaryFile = diaryFile; end

if ~isfield(globalOpt, 'saveName') || isempty(globalOpt.saveName)
    model.saveName = vargplvmWriteResult([], model.type, globalOpt.dataSetName, globalOpt.experimentNo, [], globalOpt.saveModelDir);
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



% Force kernel computations
params = svargplvmExtractParam(model);
model = svargplvmExpandParam(model, params);

if ~isempty(globalOpt.dynamicsConstrainType)
    fprintf('# Median of vardist. covars: %d \n',median(median(model.vardist.covars)));
    fprintf('# Min of vardist. covars: %d \n',min(min(model.vardist.covars)));
    fprintf('# Max of vardist. covars: %d \n',max(max(model.vardist.covars)));
end


%%

if trainModel
    model = svargplvmOptimiseModel(model);
    svargplvmShowScales(model)
end

return


%% PREDICTIONS

% obsMod = 1 means that the 1st modality is considered to be observed and
% the other modality is considered to be the unobserved one, during test
% time.
obsMod = 1; % one of the involved sub-models (the one for which we have the data)
infMod = setdiff(1:2, obsMod);

[sharedDims, privateDims] = svargplvmFindSharedDims(model);

% We can either use a totally new test set (a matrix Yts must be created
% where e.g. Yts{1} and Yts{2} are the two test modalities) or using the
% training set. In the later case, we solve the correspondance problem, ie
% for a given y we find the K most similar z's in the other modality (in
% the code, K = numberOfNN). See the Yale faces example.
fprintf('\n# PREDICTIONS: \n\n');
if ~exist('testOnTraining')
    testOnTraining=0;
end

numberOfNN = 5;

%--------------------%
svargplvmPredictions %---- Script returning: ZpredMuAll and mini(the indices for NN)
%--------------------%

