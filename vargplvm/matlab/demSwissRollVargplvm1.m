% DEMSWISSROLLVARGPLVM1 Run variational GPLVM on swiss roll data.

% VARGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'swissRoll';
experimentNo = 1;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
options = vargplvmOptions('dtcvar');
options.kern = {'rbfard2', 'bias', 'white'};
options.numActive = 100; 
%options.initX = 'ppca';
%options.scale2var1 = 1; % scale data to have variance 1
options.tieParam = 'tied';  

options.optimiser = 'scg';
latentDim = 3;
d = size(Y, 2);

% demo using the variational inference method for the gplvm model
model = vargplvmCreate(latentDim, d, Y, options);

model = vargplvmParamInit(model, model.m, model.X); 
%model.vardist.covars = 10*ones(model.N,model.q) + 0.001*randn(model.N,model.q);

% Optimise the model.
iters = 15; % Default: 1500
display = 1;

model = vargplvmOptimise(model, display, iters);

% Save the results.
modelWriteResult(model, dataSetName, experimentNo);

if exist('printDiagram') & printDiagram
   % order wrt to the inputScales
  mm = vargplvmReduceModel(model,2);
  lvmPrintPlot(mm, lbls, dataSetName, experimentNo);
end


% Load the results and display dynamically.
lvmResultsDynamic(model.type, dataSetName, experimentNo, 'plot3', model.y)


lvmScatterPlotColor(mm, model.y(:, 2)); %%% ?

% See also: plotting with scatter3

