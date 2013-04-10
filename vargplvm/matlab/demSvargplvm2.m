% DEMSVARGPLVM2 A generic demo for the MRD method.
%
% COPYRIGHT: Andreas C. Damianou, 2012
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
demSvargplvm2
%}

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

if ~exist('Yall')
    error('# This demo requires that a cell-array Yall is a priori created. This array should contain one dataset per cell.')
end
M = length(Yall);

% This script initialises the options structure 'globalOpt'.
svargplvm_init;

%-- Load datasets
for i=1:M
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
clear('Y')


for i=2:M
    if N{i} ~= N{i-1}
        error('The number of observations in each dataset must be the same!');
    end
end



%--- Create model

options = svargplvmOptions(Ytr, globalOpt);



if ~isempty(globalOpt.dynamicsConstrainType)
    for i=1:M
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


if ~isfield(globalOpt, 'saveName') || isempty(globalOpt.saveName)
    modelType = model.type;
    modelType(1) = upper(modelType(1));
    globalOpt.saveName = ['dem' modelType num2str(globalOpt.experimentNo) '.mat'];
end
model.saveName = globalOpt.saveName;
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

model = svargplvmOptimiseModel(model);

svargplvmShowScales(model)
