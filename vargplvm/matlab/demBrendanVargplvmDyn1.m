% DEMBRENDANVARGPLVMDYN1 Run variational GPLVM on Brendan face data.

% VARGPLVM
%clear;

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'brendan';

clear timeStamps; % in case it's left from a previous experiment
% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Define constants (in a manner that allows other scripts to parametrize
% this one.
if ~exist('expNo')
    expNo = 404;
end
if ~exist('itNo')
    itNo = 50; % Default: 2000
end
if ~exist('indPoints')
    indPoints = 50; % Default: 50
end
if ~exist('latentDim')
    latentDim = 30; % Default: 30
end
if ~exist('dynUsed')
    dynUsed = 1;
end
if ~exist('dataToKeep')
    dataToKeep = size(Y,1);
end
if dataToKeep == -1
    dataToKeep = size(Y,1);
end

experimentNo = expNo;
printDiagram = 1;
displayDynamic = 0;


% For fewer data
if dataToKeep < size(Y,1)
    fprintf(1,'# Using only a subset of %d datapoints...\n', dataToKeep);
    Y = Y(1:dataToKeep,:);
end

% Set up model
options = vargplvmOptions('dtcvar');
options.kern = {'rbfard2', 'white'};
options.numActive = indPoints; % Default: 100
options.scale2var1 = 1; % scale data to have variance 1
%options.tieParam = 'tied';  

options.optimiser = 'scg';
%latentDim = 10; % Default: 30
d = size(Y, 2);

% create the model
fprintf(1,'# Creating the model...\n');
model = vargplvmCreate(latentDim, d, Y, options);
%
model = vargplvmParamInit(model, model.m, model.X); 


%-------- Add dynamics to the model -----
if dynUsed
    model = addDefaultVargpTimeDynamics(model);
end
%------------------ 

if ~isempty(model.dynamics)
% A better way to initialize the  kernel hyperparameter,
% especially lengthscale, should be to learn it by running few iterations of
% GP regression marginal likelihood maximization given as data the PCA output
% (the default value is jsut reasonable and it sets the inverse lenthscale to quite 
% small value so the GP dynamic prior is weaker (less smoothed) at
% initializations
model.dynamics.kern.comp{1}.inverseWidth = 20./(((max(model.dynamics.t)-min(model.dynamics.t))).^2);
params = vargplvmExtractParam(model);
model = vargplvmExpandParam(model, params);    
% Initialize barmu     
initFunc = str2func([options.initX 'Embed']);
X = initFunc(model.m, model.q);
vX = var(X); 
for q=1:model.q
    Lkt = chol(model.dynamics.Kt + 0.01*vX(q)*eye(model.N))';    
    % barmu = inv(Kt + s*I)*X, so that  mu = Kt*barmu =  Kt*inv(Kt +
    % s*I)*X, which is jsut the GP prediction, is temporally smoothed version 
    % of the PCA latent variables X (for the data Y)
    model.dynamics.vardist.means(:,q) = Lkt'\(Lkt\X(:,q));    
end
params = vargplvmExtractParam(model);
model = vargplvmExpandParam(model, params);
% inducing point need to initilize based on model.vardist.means
perm = randperm(model.k); 
model.X_u = model.vardist.means(perm(1:model.k),:);
params = vargplvmExtractParam(model);
model = vargplvmExpandParam(model, params);
end
modelInit = model;


% Optimise the model.
iters = itNo; % DEfault: 1500
display = 1;

fprintf(1,'# Optimising the model...\n');
model.learnBeta = 1;
model = vargplvmOptimise(model, display, iters);

% Save the results.
fprintf(1,'# Saving the model...\n');
modelWriteResult(model, dataSetName, experimentNo);

if exist('printDiagram') & printDiagram
  lvmPrintPlot(model, lbls, dataSetName, experimentNo);
end

if exist('displayDynamic') & displayDynamic
    % load connectivity matrix
    %[void, connect] = mocapLoadTextData('run1');
    % Load the results and display dynamically.
    lvmResultsDynamic(model.type, dataSetName, experimentNo, 'image', [20 28],1, 0, 1)
end
