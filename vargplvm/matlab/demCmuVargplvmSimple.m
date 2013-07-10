% DEMCMU35GPLVMVARGPLVMSIMPLE Run variational GPLVM with dynamics on a single
% sequence of the CMU35 data.

% VARGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

% Define constants
experimentNo = 1000;
latentDim = 9; % this is Q, the number of latent dimensions
indPoints = 100; % number of inducing points
% dynamicKern = {'matern32', 'bias', 'white'};
dynamicKern = {'rbf', 'bias', 'white'}; % kernel k_t for the GP for x(t)
initX ='ppca'; % initialize latent space with ppca

% load data
dataSetName = 'cmu35gplvm';
[Ytmp, lbls, Y, lblstest] = lvmLoadData(dataSetName);

% Set training and test sets
indicesTraining = [1:40 61:size(Y,1)];
indicesTest = setdiff(1:size(Y,1),indicesTraining);
Ytr = Y(indicesTraining,:);
Ytest = Y(indicesTest,:);

% Corresponding timestamps (artificial and equally spaced for this demo)
t = linspace(0, 2*pi, size(Y, 1)+1)';
t = t(1:end-1, 1);
timeStampsTraining = t(indicesTraining);
timeStampsTest = t(indicesTest);

Y = Ytr; clear('Ytr','Ytmp');


% Set up model
options = vargplvmOptions('dtcvar');
options.kern = {'rbfard2', 'bias', 'white'}; % Kernel k_x for the GP for f(x)
options.numActive = indPoints;
options.optimiser = 'scg';

d = size(Y, 2);
fprintf(1,'# Creating the model...\n');
model = vargplvmCreate(latentDim, d, Y, options);
model = vargplvmParamInit(model, model.m, model.X);
model.beta=1/(0.01*var(model.m(:)));
model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));

%-------- Add dynamics to the model -----
optionsDyn.type = 'vargpTime';
optionsDyn.t=timeStampsTraining;
optionsDyn.inverseWidth=30;
optionsDyn.initX = initX;

% Dynamic kernel:
kern = kernCreate(t, dynamicKern);
% The following is related to the expected number of
% zero-crossings.(larger inv.width numerator, rougher func)
if ~strcmp(kern.comp{1}.type,'ou')
    kern.comp{1}.inverseWidth = optionsDyn.inverseWidth./(((max(t)-min(t))).^2);
    kern.comp{1}.variance = 1;
end
optionsDyn.kern = kern;

% Fill in with default values whatever is not already set
optionsDyn = vargplvmOptionsDyn(optionsDyn);
model = vargplvmAddDynamics(model, 'vargpTime', optionsDyn, optionsDyn.t, 0, 0,optionsDyn.seq);

fprintf(1,'# Further calibration of the initial parameters...\n');
model = vargplvmInitDynamics(model,optionsDyn);
model.vardist.parallel=1;
% do not learn beta for few iterations for intitilization
model.learnBeta = 0;
display = 1;
fprintf(1,'# Intitiliazing the model (fixed beta) %d iterations...\n',100);
model = vargplvmOptimise(model, display, 100);
disp('# Saving model after optimising beta...')
modelWriteResult(model, dataSetName, experimentNo);

% Optimise the model.
model.learnBeta = 1;
iters = 1000; % Default: 1000
fprintf(1,'# Optimising the model for %d iterations...\n',iters);
model = vargplvmOptimise(model, display, iters);
% Save the results.
fprintf(1,'# Saving model after doing %d iterations\n',iters)
modelWriteResult(model, dataSetName, experimentNo);

% See the final lengthscales (see how some dimensions are switched-off).
 bar(model.kern.comp{1}.inputScales)

%% ----------------- PREDICTIONS -----------------------------
fprintf('# Only times prediction...\n');
% Prediction using the only information in the test time points
[Testmeans2 Testcovars2] = vargplvmPredictPoint(model.dynamics, timeStampsTest);
Varmu2 = vargplvmPosteriorMeanVar(model, Testmeans2, Testcovars2);

errorPredictions = mean(abs(Varmu2(:) - Ytest(:)));
fprintf(1,'# Predictions Error:%d\n', errorPredictions);


plot(Varmu2(1,:),'r'),hold on, plot(Ytest(1,:),'g')

skel = acclaimReadSkel('35.asf');
[tmpchan, skel] = acclaimLoadChannels('35_01.amc', skel);
channels{1} = demCmu35VargplvmLoadChannels(Ytest,skel);
channels{2} = demCmu35VargplvmLoadChannels(Varmu2,skel);
skelPlayData2(skel, channels,1/15,{'Ytest','Varmu2'});
