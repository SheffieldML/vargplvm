% DEMOVARGPLVM1 Description ...

% VARGPLVM



% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'oil';
experimentNo = 1;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
options = vargplvmOptions('dtcvar');
options.kern = {'rbfard2'};
options.numActive = 10; 

options.optimiser = 'scg';
latentDim = 2;
d = size(Y, 2);

%  demo using the variational inferecne method for the gplvm model 
model = vargplvmCreate(latentDim, d, Y, options);

params = vargplvmExtractParam(model);
model2 = vargplvmExpandParam(model,params);
