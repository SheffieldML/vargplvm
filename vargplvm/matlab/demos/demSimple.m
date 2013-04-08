
% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);


% Constants:
%dataSetName = 'brendan';
%dataSetName = 'eth-car1';
dataSetName ='eth-tomato1';
experimentNo = 404;
itNo = 200;
initVardistIters = 200;


[Y tmp] = lvmLoadData(dataSetName);
height = tmp(1);
width = tmp(2);
clear tmp;

dims = size(Y,2);


% Split into training and test set.
Ytr = Y;

%Ytr = Ytr(:,1); dims=1;%%%%

% Set up model
options = vargplvmOptions('dtcvar');
options.kern = {'rbfard2', 'bias', 'white'}; % mapping kernel
%options.kern = 'rbfardjit'; % mapping kernel
options.numActive = 30; % ind. points
options.optimiser = 'scg';
latentDim = 10; % Q %%%% DEF: 10
d = size(Ytr, 2);


% demo using the variational inference method for the gplvm model
fprintf(1,'# Creating the model...\n');
% scale = std(Ytr);
% scale(find(scale==0)) = 1;
%options.scaleVal = mean(std(Ytr));
options.scaleVal = sqrt(var(Ytr(:)));

model = vargplvmCreate(latentDim, d, Ytr, options);



% Temporary: in this demo there should always exist the mOrig field
if ~isfield(model, 'mOrig')
    model.mOrig = model.m;
end

% Initialise the model
model = vargplvmParamInit(model, model.mOrig, model.X);
invWidthMult = 5;
model.kern.comp{1}.inputScales = invWidthMult./(((max(model.X)-min(model.X))).^2); % Default 5
params = vargplvmExtractParam(model);
model = vargplvmExpandParam(model, params);

model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));
model.beta=1/(0.01*var(model.mOrig(:)));

modelInit = model;


%-
capName = dataSetName;
capName(1) = upper(capName(1));
modelType = model.type;
modelType(1) = upper(modelType(1));
fileToSave = ['dem' capName modelType num2str(experimentNo) '.mat'];
%

display = 1;

% Initialise the var. distr
if initVardistIters ~=0
    model.initVardist = 1;
    fprintf(1,'# Intitiliazing the model (fixed beta)...\n');
    model = vargplvmOptimise(model, display, initVardistIters); % Default: 20
    fprintf(1,'1/b = %.4d\n',1/model.beta);
    
    model.initVardist = 0;
end

model.learnBeta = 1;
model.iters = 0;

prunedModelInit = vargplvmPruneModel(modelInit);
clear modelInit


% Optimise the model.
for i=1:length(itNo)
    iters = itNo(i); % default: 2000
    fprintf(1,'\n# Optimising the model for %d iterations (session %d)...\n',iters,i);
    model = vargplvmOptimise(model, display, iters);
    model.iters = model.iters + iters;
    
    fprintf(1,'# 1/b=%.4f\n var(m)=%.4f\n',1/model.beta, var(model.mOrig(:)));
    
    % Save model
    fprintf(1,'# Saving %s\n',fileToSave);
    prunedModel = vargplvmPruneModel(model);
    % prunedModelTr = vargplvmPruneModel(modelTr);
    save(fileToSave, 'prunedModel', 'prunedModelInit');
end
prunedModelTr = prunedModel;
save(fileToSave, 'prunedModel', 'prunedModelInit', 'prunedModelTr');

model.m = model.mOrig;
[modelP, newHeight, newWidth] = vargplvmReduceVidModel(model, height, width, 4,4);
lvmVisualise(modelP, [], 'imageVisualise', 'imageModify', [newHeight newWidth],0,0,1);
clear modelP


