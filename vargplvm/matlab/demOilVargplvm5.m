% DEMOILVARGPLVM5 Run variational GPLVM on oil data.

% VARGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'oil';
experimentNo = 5;
printDiagram = 1;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
options = vargplvmOptions('dtcvar');
%options.kern = {'linard2', 'bias', 'white'};
options.kern = {'rbfard2','bias','white'};
options.numActive = 50; 
%options.tieParam = 'tied';  

options.optimiser = 'scg';
latentDim = 10;
d = size(Y, 2);

% demo using the variational inference method for the gplvm model
model = vargplvmCreate(latentDim, d, Y, options);
%
model = vargplvmParamInit(model, model.m, model.X); 
model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));
model.learnBeta=1;
model.beta = 1/((1/100 * var(model.m(:))));
% Optimise the model.
iters = 2500;
display = 1;

model = vargplvmOptimise(model, display, iters);
model.info = info;
capName = dataSetName;
capName(1) = upper(capName(1));
modelType = model.type;
modelType(1) = upper(modelType(1));
save(['dem' capName modelType num2str(experimentNo) '.mat'], 'model');

% order wrt to the inputScales 
mm = vargplvmReduceModel(model,2);
%% plot the two largest twe latent dimensions 
if exist('printDiagram') & printDiagram
  lvmPrintPlot(mm, lbls, capName, experimentNo);
end
