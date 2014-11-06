% DEMBRENDANVARGPLVM2 Run variational GPLVM on Brendan face data.

% VARGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'brendan';
experimentNo = 2;
printDiagram = 1;

% load data
[Y, lbls] = lvmLoadData('brendan');

% Set up model
options = vargplvmOptions('dtcvar');
options.kern = {'rbfard2', 'white'};
options.numActive = 50; 
%options.scale2var1 = 0; % scale data to have variance 1
%options.tieParam = 'tied';  

options.optimiser = 'scg';
latentDim = 30;
d = size(Y, 2);

% demo using the variational inference method for the gplvm model
model = vargplvmCreate(latentDim, d, Y, options);
%
model = vargplvmParamInit(model, model.m, model.X); 

% Optimise the model.
iters = 1500;
display = 1;

model = vargplvmOptimise(model, display, iters);


capName = dataSetName;
capName(1) = upper(capName(1));
modelType = model.type;
modelType(1) = upper(modelType(1));
save(['dem' capName modelType num2str(experimentNo) '.mat'], 'model');

% order wrt to the inputScales 
model = vargplvmOrderLatentDims(model);

% plot the two largest latent dimensions 
if exist('printDiagram') & printDiagram
  lvmPrintPlot(model, lbls, capName, experimentNo);
end

%% Load the results and display dynamically.
%fgplvmResultsDynamic(dataSetName, experimentNo, 'image', [20 28], 1, 0, 1)

%% compute the nearest neighbours errors in latent space.
%errors = fgplvmNearestNeighbour(model, lbls);
