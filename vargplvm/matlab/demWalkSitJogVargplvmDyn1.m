%
%	DEMWALKSITJOGVARGPLVMDYN1: Demo of for the walksitjog mocap dataset testing the
% dynamics of vargplvm.
%	
  
% VARGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'walkSitJog';

clear timeStamps; % in case it's left from a previous experiments
% load data; the function called here is like lvmLoadData but also loads
% timeStamps.
[Y timeStamps] = loadWalkSitJog;
lbls = 'connect';
%%%

if ~exist('expNo')
    expNo = 404;
end
if ~exist('itNo')
    itNo = 40; % Default: 2000
end
if ~exist('indPoints')
    indPoints = 10; % Default: 50
end
if ~exist('latentDim')
    latentDim = 10;
end
if ~exist('dynUsed')
    dynUsed = 1;
end
if ~exist('dataToKeep')
    dataToKeep = 1000
end
if dataToKeep == -1
    dataToKeep = size(Y,1);
end

experimentNo = expNo;
printDiagram = 1;
displayDynamic = 0;


% for fewer data
if dataToKeep < size(Y,1)
    fprintf(1,'# Using only a subset of %d datapoints...\n', dataToKeep);
    Y = Y(1:dataToKeep,:);
    timeStamps = timeStamps(1:dataToKeep,:);
end


% Set up model
options = vargplvmOptions('dtcvar');
options.kern = {'rbfard2', 'bias', 'white'};
options.numActive = indPoints; % Default: 50
%options.tieParam = 'tied';  

options.optimiser = 'scg';
%options.optimiser = 'optimiMinimize';
%latentDim = 10; % Default: 10
d = size(Y, 2);

% demo using the variational inference method for the gplvm model
fprintf(1,'# Creating the model...\n');
model = vargplvmCreate(latentDim, d, Y, options);
%
model = vargplvmParamInit(model, model.m, model.X); 
model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));



%-------- Add dynamics to the model -----
if dynUsed
    model = addDefaultVargpTimeDynamics(model, timeStamps);
end
%------------------ 


% Optimise the model.
iters = itNo; % Default: 2000
display = 1;
fprintf('# Optimising the model...\n');
model = vargplvmOptimise(model, display, iters);

fprintf(1,'# Saving the model...\n');
capName = dataSetName;;
capName(1) = upper(capName(1));
modelType = model.type;
modelType(1) = upper(modelType(1));
save(['dem' capName modelType num2str(experimentNo) '.mat'], 'model');

% order wrt to the inputScales 
mm = vargplvmReduceModel(model,2);
% plot the two largest twe latent dimensions 
if exist('printDiagram') & printDiagram
  lvmPrintPlot(mm, lbls, capName, experimentNo);
end
%lvmResultsDynamic(model.type, dataSetName, experimentNo, 'robotWireless',
%'vector')
