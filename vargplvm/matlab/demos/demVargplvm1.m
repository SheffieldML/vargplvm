% Examples:
% dataSetName = 'hierarchical/demHighFiveHgplvm1'; dataSetField = 'YA';
%      demVargplvm1

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);


vargplvm_init;

%--- Load data
% There should exist either a dataset Y already loaded, or a dataSetName
% with the name of the dataset.
if ~isempty(globalOpt.dataSetName)
    try
        Y = vargplvmLoadData(globalOpt.dataSetName,[],[],globalOpt.dataSetField);
    catch
        error('Could not load dataset')
    end
elseif ~exist('Y')
    error('dataSetName or Y should be defined.')
end


%%

% Set up model
options = vargplvmOptions('dtcvar');
options.kern = globalOpt.mappingKern;
options.numActive = globalOpt.indPoints;
if options.numActive == -1
    options.numActive = size(Y,1);
end
options.optimiser = 'scg2';
options.initX = globalOpt.initX;
options.enableDgtN = globalOpt.DgtN;
%options.scaleVal = sqrt(var(Y(:)));
if globalOpt.fixInd
    options.fixInducing=1;
    options.fixIndices=1:size(Y,1);
end
options.scale2var1 = globalOpt.scale2var1;
fprintf(1,'# Creating the model...\n');
model = vargplvmCreate(globalOpt.latentDim, size(Y,2), Y, options);
model = vargplvmModelInit(model, globalOpt);
model.globalOpt = globalOpt;

modelInit = model; %%%


%%
model = vargplvmOptimiseModel(model, true, true, {globalOpt.initVardistIters,globalOpt.itNo});

if ~isempty(globalOpt.saveName)
    save(globalOpt.saveName, model);
end
