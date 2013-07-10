% DEMOIL100POINT Description

% VARGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'oil100';
experimentNo = 100;
printDiagram = 1;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
options = vargplvmOptions('dtcvar');
options.kern = {'linard2', 'bias', 'white'};
options.numActive = 50; 
%options.tieParam = 'tied';  

options.optimiser = 'scg';
latentDim = 12;
d = size(Y, 2);

% demo using the variational inference method for the gplvm model
model = vargplvmCreate(latentDim, d, Y, options);
%
model = vargplvmParamInit(model, model.m, model.X); 

iters = 100;
display = 1;

model = vargplvmOptimise(model, display, iters);

ll = vargplvmLogLikelihood(model);
g = vargplvmLogLikeGradients(model);
gmeans = reshape(g(1:(model.vardist.numData*model.vardist.latentDimension)), model.vardist.numData, model.vardist.latentDimension);
gcovars = reshape(g((model.vardist.numData*model.vardist.latentDimension+1):2*model.vardist.numData*model.vardist.latentDimension), model.vardist.numData, model.vardist.latentDimension);

model1 = model;
model1.m = model.m(2:end,:)
model1.vardist.means = model1.vardist.means(2:end,:); 
model1.vardist.covars = model1.vardist.covars(2:end,:);
model1.X = model.X(2:end,:); 
model1.N = model1.N-1;
model1.vardist.numData = model1.vardist.numData-1;
model1.Psi1 =  kernVardistPsi1Compute(model1.kern, model1.vardist, model1.X_u);
model1.Psi2 =  kernVardistPsi2Compute(model1.kern, model1.vardist, model1.X_u);
model1.Psi0 =  kernVardistPsi0Compute(model1.kern, model1.vardist);

vardistx = vardistCreate(model.vardist.means(1,:), model.q, 'gaussian');
vardistx.means = model.vardist.means(1,:); 
vardistx.covars = model.vardist.covars(1,:);
indexPresent = 1:model.q;
ll1 = vargplvmPointLogLikelihood(model1, vardistx, model.y(1,indexPresent), indexPresent);
g1 = vargplvmPointLogLikeGradient(model1, vardistx, model.y(1,indexPresent), indexPresent);


if 0 
% Optimise the model.
iters = 2000;
display = 1;

model = vargplvmOptimise(model, display, iters);

capName = dataSetName;;
capName(1) = upper(capName(1));
modelType = model.type;
modelType(1) = upper(modelType(1));
save(['dem' capName modelType num2str(experimentNo) '.mat'], 'model');

% order wrt to the inputScales 
model = vargplvmOrderLatentDims(model);
% plot the two largest latent dimensions 
model.X = model.vardist.means(:,[2 1]); % for the oil plot the largest in the vertical axis
if exist('printDiagram') & printDiagram
  lvmPrintPlot(model, lbls, capName, experimentNo);
end
end
