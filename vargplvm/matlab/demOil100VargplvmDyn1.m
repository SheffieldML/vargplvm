% DEMOILVARGPLVM1 Run variational GPLVM on oil data.

% VARGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

if ~exist('vardistCovarsMult'), vardistCovarsMult = 0.5; end

dataSetName = 'oil100';
if ~exist('experimentNo'), experimentNo = 5; end
printDiagram = 1;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
options = vargplvmOptions('dtcvar');
%options.kern = {'rbfard2', 'bias', 'white'};
options.kern = 'rbfardjit';
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



   if dynUsed
        timeStampsTraining = [1:model.N]';
        t=timeStampsTraining;
        optionsDyn.inverseWidth=100;
        fprintf(1,'# Adding dynamics to the model...\n');
        optionsDyn.type = 'vargpTime';
        optionsDyn.t=timeStampsTraining;
        %optionsDyn.inverseWidth=invWidthMultDyn; % Default: 100
        %optionsDyn.testReoptimise = testReoptimise;
        dynamicKern = {'rbf','white','bias'};        
        kern = kernCreate(t, dynamicKern); % Default: {'rbf','white','bias'}
        
        
        %-- Initialize each element of the compound kernel (optional but
        % recommended)
        % ATTENTION: For the gradients we assume that the base kernel (rbf,
        % matern etc) must be the FIRST one and if a second base kernel
        % (not white or bias) exist must be the LAST one!!!!!!!!!!!!!!
        vargplvmInitDynKernel;
        
        optionsDyn.learnVariance = 1; %% NEW
        optionsDyn.kern = kern;
        optionsDyn.vardistCovars = vardistCovarsMult; % 0.23 gives true vardist.covars around 0.5 (DEFAULT: 0.23) for the ocean dataset
        
        % Fill in with default values whatever is not already set
        optionsDyn = vargplvmOptionsDyn(optionsDyn);
        model = vargplvmAddDynamics(model, 'vargpTime', optionsDyn, optionsDyn.t, 0, 0,optionsDyn.seq);
        
        fprintf(1,'# Further calibration of the initial values...\n');
        model = vargplvmInitDynamics(model,optionsDyn);
   end

% Optimise the model.
iters = 1800;
display = 1;

model.beta = 1/((1/100 * var(model.m(:))));
model.learnBeta = false; model.learnSigmaf = false; model.initVardist = true;
model = vargplvmOptimise(model, display, 500);
model.learnBeta = true; model.learnSigmaf = true; model.initVardist = false;


model = vargplvmOptimise(model, display, iters);

capName = dataSetName;
capName(1) = upper(capName(1));
modelType = model.type;
modelType(1) = upper(modelType(1));
save(['dem' capName modelType num2str(experimentNo) '.mat'], 'model');

% order wrt to the inputScales 
mm = vargplvmReduceModel2(model,2);
%% plot the two largest twe latent dimensions 
if exist('printDiagram') & printDiagram
    try
  lvmPrintPlot(mm, lbls, capName, experimentNo);
    catch e
    end
  figure
  bar(model.kern.inputScales);
end
errors = fgplvmNearestNeighbour(mm, lbls);
fprintf('# Vargplvm errors in the 2-D projection: %d\n', errors)

%%
mmm = vargplvmReduceModel2(model, length(vargplvmRetainedScales(model)));
fprintf('# Bound in the reduced model: %.4f\n', -vargplvmLogLikelihood(mmm));