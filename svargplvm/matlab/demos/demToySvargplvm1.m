% DEMTOYSVARGPLVM1 A simple demo of Manifold Relevance Determination.
% DESC This is a simple demo which uses a toy dataset that considers two
% modalities: each of the modalities has one private 1-D signal and they
% both share a 1-D shared signal. The final dataset is mapped to 10-D. 
% The demo allows training a MRD model so that the learned weights segment
% the output space and recover the true signals. 
% 
% The demo demonstrates a basic options (many more are generally
% available):
% experiment with different latent space initialisations, different kernels
% (linear, non-linear, etc) and with and without constraining the latent
% space with dynamics. Visualise the learned spaces and scales, see the
% difference in smoothness when dynamics are used/not used, compare the
% variational bounds, etc. 
% 
% Examples of running the demo:
% a) Non-dynamical case, initialise X separately for each modality and then concatenate:
% --------------
% >> latentDimPerModel = 4;
% >> initial_X = 'separately';
% >> demToySvargplvm1
%
% b) Non-dynamical case, initialise X once for the outputs that are first
% concatenated
% --------------
% >> latentDim = 6;
% >> initial_X = 'concatenated'; 
% >> demToySvargplvm1
%
% c) Dynamics: 
% The above runs can be combined with dynamics (requires more iterations to converge): 
% >> initVardistIters = 600; 
% >> itNo = 1000;
% >> dynamicKern = {'rbf', 'white', 'bias'}; % SEE KERN toolbox for kernels
% >> dynamicsConstrainType = {'time'};
% >> demToySvargplvm1
%
%
% For more options / variations, check/alter the variables in the first section 
% the demo as well as the initialisation function svargplvm_init.m which
% lists all possible options.
%
% COPYRIGHT: Andreas C. Damianou, 2012
%
% SHEFFIELDML

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

%% ---  High level options for the demo
%-- Mapping kernel: One cell for each modality. Compound kernels (cells
%themselves) allowed, e.g. {k1, k2} means k1+k2.
if ~exist('baseKern', 'var'), baseKern = {{'linard2','white', 'bias'},{'linard2','white', 'bias'}}; end
%-- # iterations of initialisation of the variational distribution
if ~exist('initVardistIters', 'var'), initVardistIters = 200; end
%-- # iterations of optimisation
if ~exist('itNo', 'var'), itNo = 200; end
%-- Initialisation of latent space:
if ~exist('latentDimPerModel', 'var'), latentDimPerModel = 4; end % Only used if initial_X is 'separately'
if ~exist('latentDim', 'var'), latentDim = 6; end % Only used if initial_X is 'concatenated'
if ~exist('initial_X', 'var'), initial_X = 'separately';  end


%% Create toy dataset
alpha = linspace(0,4*pi,100);
% Scale and center data
Z1 = scaleData(cos(alpha)', 2);
Z2 = scaleData(sin(alpha)', 2);
Z3 = scaleData((cos(alpha)').^2, 2); % OR: 2*cos(2*alpha)' + 2*sin(2*alpha)'
noiseLevel = 0.1; % Default: 0.1
% Map 1-D to 10-D and add some noise
Z2p = Z2*rand(1,10);
Z2p = Z2p + noiseLevel.*randn(size(Z2p));
Z1p = Z1*rand(1,10);
Z1p = Z1p + noiseLevel.*randn(size(Z1p));
Z3p = Z3*rand(1,10);%
Z3p = Z3p + noiseLevel.*randn(size(Z3p));%
%---
numSharedDims = 5;
Z1p(:,1:numSharedDims) = Z3p(:,1:numSharedDims);
Z2p(:,1:numSharedDims) = Z3p(:,1:numSharedDims);
bar(pca([Z1p Z2p]))
title(sprintf('PCA scales in the final %d-dimensional dataset.', size(Z1p,2)+size(Z2p,2)))
%---
Yall{1} = Z1p;
Yall{2} = Z2p;
dataSetNames={'fols_cos', 'fols_sin'};
M = 2; % number of modalities

% This script initialises the options structure 'globalOpt'.
svargplvm_init;

%% Split data into training and test sets.
for i=1:M
    Y = Yall{i};
    dims{i} = size(Y,2);
    N{i} = size(Y,1);
    indTr = globalOpt.indTr;
    if indTr == -1, indTr = 1:N{i}; end
    if ~exist('Yts', 'var')
        indTs = setdiff(1:size(Y,1), indTr);
        Yts{i} = Y(indTs,:);
    end
    Ytr{i} = Y(indTr,:);
end
t = linspace(0, 2*pi, size(Y, 1)+1)'; t = t(1:end-1, 1);
timeStampsTraining = t(indTr,1); %timeStampsTest = t(indTs,1);
clear('Y')
for i=2:M
    if N{i} ~= N{i-1}
        error('The number of observations in each dataset must be the same!');
    end
end

%% -- Create model
options = svargplvmOptions(Ytr, globalOpt);

% Allow for constraining the latent space with a prior which couples it.
% This can e.g be a termporal one (type: 'vargpTime') for which timestamps can
% be given, or a labelset (type: 'labelsTrain'), for which the labels must
% be given
if ~isempty(globalOpt.dynamicsConstrainType)
    optionsDyn.type = 'vargpTime';
    optionsDyn.inverseWidth=30;
    optionsDyn.initX = globalOpt.initX;
    optionsDyn.constrainType = globalOpt.dynamicsConstrainType;
    if exist('timeStampsTraining', 'var')
        optionsDyn.t = timeStampsTraining;
    end
    if exist('labelsTrain', 'var') && ~isempty(labelsTrain)
        optionsDyn.labels = labelsTrain;
    end
else
    optionsDyn= [];
end

model = svargplvmModelCreate(Ytr, globalOpt, options, optionsDyn);

if ~isfield(globalOpt, 'saveName') || isempty(globalOpt.saveName)
    modelType = model.type;
    modelType(1) = upper(modelType(1));
    globalOpt.saveName = ['dem' modelType num2str(globalOpt.experimentNo) '.mat'];
end
model.saveName = globalOpt.saveName;
model.globalOpt = globalOpt;
model.options = options;

%-- Define what level of parallelism to use (w.r.t submodels or/and w.r.t
% datapoints).
%{
fprintf('# Parallel computations w.r.t the submodels!\n');
model.parallel = 1;
model = svargplvmPropagateField(model,'parallel', 1);
%
fprintf('# Parallel computations w.r.t the datapoints!\n');
model.vardist.parallel = 1;
for i=1:model.numModels
    model.comp{i}.vardist.parallel = 1;
end
%}

% Force kernel computations
params = svargplvmExtractParam(model);
model = svargplvmExpandParam(model, params);

%% Optimise model
model = svargplvmOptimiseModel(model);

%% Results
fprintf('# Finished!\n')
fprintf('# Variational lower bound: %.6f\n', svargplvmLogLikelihood(model))
figure
svargplvmShowScales(model); title('Learned scales for the two models')
[sharedDims, privateDims] = svargplvmFindSharedDims(model,0.05,true);
%
figure
subplot(1,2,1)
plot(model.vardist.means(:,sharedDims), 'x-');
title('Recovered signal for shared dimensions.')
subplot(1,2,2)
plot(Z3, 'x-'); title('True (and noiseless) shared signal.');
%
for i = 1:length(privateDims)
figure
subplot(1,2,1)
plot(model.vardist.means(:,privateDims{i}), 'x-');
title(sprintf('Recovered signal for private dimensions of model %d.', i))
subplot(1,2,2)
if i == 1
    plot(Z1, 'x-'); 
else
    plot(Z2, 'x-'); 
end
title(sprintf('True (and noiseless) private signal for modality %d.', i));
end

if ~isempty(globalOpt.dynamicsConstrainType)
    figure
    imagesc(model.dynamics.Kt); title('Latent space prior cov. matrix.')
end