function defaults = defaults_svargplvm
% COPYRIGHT: Andreas C. Damianou, 2012, 2013
% VARGPLVM

% Fields Stored as: {fieldName1, value1, fieldName2, value2, ... }
%

defaults = {};
%% Default values for the model
defaults{end+1} = 'experimentNo'; defaults{end+1} = 404;
% Set this value to 'noSave' if you don't want to save the model
defaults{end+1} = 'saveName'; defaults{end+1} = [];
defaults{end+1} = 'itNo'; defaults{end+1} = [500 2000];
defaults{end+1} = 'indPoints'; defaults{end+1} = 100;
defaults{end+1} = 'latentDim'; defaults{end+1} = 10;
% This should be a cell array of length M (num. of modalities) if the
% field initial_X is set to 'custom'
defaults{end+1} = 'latentDimPerModel'; defaults{end+1} = 3;
defaults{end+1} = 'initVardistIters'; defaults{end+1} = 300;
% How to initialise X:
%   'separately', means apply the initX function to each
% of the datasets and then concatenate.
%   'concatenated', means first concatenate
% the datasets and then apply the 'initX' function.
%   'custom' is like the "separately" options, but it implies that
% latentDimPerModel is a cell specifying how many dimensions to use for each submodel.
defaults{end+1} = 'initial_X'; defaults{end+1} = 'separately'; % other options: 'concatenated' or 'custom'

% How to initialise the latent space in level h based on the data of
% level h-1. This can either be a signle entry (replicated in all
% levels) or it can be a cell array of length H, i.e. different
% initialisation for each level.
% The values of this field are either strings, so that [@initX 'Embed']
% is called, or it can be a matrix of an a priori computed latent space.
% Other options for embedding:
% 'pca','isomap2','vargplvm','fgplvm','outputs'. The last requires that
% Q==D and initialises X=Y where the data is scaled so that it's 0 mean
% 1 variance
defaults{end+1} = 'initX'; defaults{end+1} = 'ppca';
defaults{end+1} = 'inputScales'; defaults{end+1} = [];
% Set to 1 to tie the inducing points with the latent vars. X
defaults{end+1} = 'fixInd'; defaults{end+1} = 0;
defaults{end+1} = 'baseKern'; defaults{end+1} = {'rbfardjit', 'rbfardjit'};   %{'rbfard2', 'white'};
defaults{end+1} = 'dynamicKern'; defaults{end+1} = {'rbf','white', 'bias'};
% Like dataToKeep, but the rest indices go to the test set.
defaults{end+1} = 'indTr'; defaults{end+1} = -1;
%%%defaults{end+1} = 'backConstraints'; defaults{end+1} = 1; % Substituted with dynamicsConstrainType
defaults{end+1} = 'vardistCovarsMult'; defaults{end+1} = 2;
% The "ideal" initial variational covariances for the dynamics layer. Since
% we do not have direct access to these covariances, we can only set the
% reparametrised covariances (lambda) which result in different var.
% covars. That's why we aim for a reasonable median.
defaults{end+1} = 'dynamicsInitCovarMedian'; defaults{end+1} = 0.18;
% The lowest and highest acceptable bounds for the above.
defaults{end+1} = 'initCovarMedianLowest'; defaults{end+1} = [];
defaults{end+1} = 'initCovarMedianHighest'; defaults{end+1} = [];
defaults{end+1} = 'dataSetNames'; defaults{end+1} = {};
defaults{end+1} = 'dataSetName'; defaults{end+1} = 'unknown';
defaults{end+1} = 'mappingInitialisation'; defaults{end+1} = 0;
defaults{end+1} = 'scale2var1'; defaults{end+1} = 0;
% Set to -1 to use all data. Set to scalar d, to only take d points from
% each class. Set to a vector D with length equal to the number of classes C,
% to take D(c) points from class c.
defaults{end+1} = 'dataPerClass'; defaults{end+1} = -1;
% Set to -1 to keep all the training data, set to a number N to only keep N
% datapoints.
defaults{end+1} = 'dataToKeep'; defaults{end+1} = -1;
% Signal to noise ratio (initialisation for model.beta).
defaults{end+1} = 'initSNR'; defaults{end+1} = 100;
% How many iterations to do to initialise the model with a static B. gplvm.
% -1 means that no such initialisation will be done.
defaults{end+1} = 'initWithStatic'; defaults{end+1} = -1;
% If initWithStatic ~= -1, this says how many iters with fixed
% beta/sigmaf to perform.
defaults{end+1} = 'initWithStaticInitVardist'; defaults{end+1} = 300;
% If initWithStatic ~= -1, this says what the initial SNR will be for the
% initial static model.
defaults{end+1} = 'initStaticSNR'; defaults{end+1} = 25;
% If true, then if also initWithStatic ~=1, we initialise the model.beta
% and model.kern based on the initialised static model.
defaults{end+1} = 'initWithStaticAll'; defaults{end+1} = false;
% Leave empty {} for no dynamics.
% Other options are: 'time', 'labels', 'data'
defaults{end+1} = 'dynamicsConstrainType'; defaults{end+1} = {};
% If set to fals, then the dynamics kernel is not initialised with
% bc_initdynKernel
defaults{end+1} = 'initDynKernel'; defaults{end+1} = 1;
% A second (probably better) way to initialise the model
defaults{end+1} = 'initDynKernel2'; defaults{end+1} = 0;
% See bc_backConstraintsModelCreate and bc_restorePrunedModel
defaults{end+1} = 'labelMatrixPower'; defaults{end+1} = 0.5;
% if not empty, the corresponding gradients of the kernel will be zero,
% i.e. not learning these elements (typically for the variance of the
% rbf/matern etc of a invcmpnd)
defaults{end+1} = 'fixedKernVarianceIndices'; defaults{end+1} = [];
% If discrKernel is 'ones', then the discriminative kernel is
% simply build based on a matrix with ones and minus ones, otherwise it is
% based on a measure on the distance of each label from the mean of each
% class.
defaults{end+1} = 'discrKernel'; defaults{end+1} = 'ones';
% Default variance for a fixedwhite kernel
defaults{end+1} = 'fixedwhiteVar'; defaults{end+1} = 1e-5;
% If set to some value, call it x, then after learning a constrained model,
% (and if the function runStaticModel is called), a static model will be
% initialised with the constrained model's latent space and learned for x iterations.
defaults{end+1} = 'runStaticModel'; defaults{end+1} = -1;
defaults{end+1} = 'runStaticModelInitIters'; defaults{end+1} = [];
% Possible values: (none, one or both of them): 'labelsYinputs' and
% 'labelsYoutputs'. If the first is there, then the Y of p(X|Y) is
% augmented with the labels as C extra dimensions, where C is the total
% number of classes (as we use 1-of-K encoding). Similarly with labelsYoutputs.
defaults{end+1} = 'dataConstraints'; defaults{end+1} = {};
% Option to only run a static model.
defaults{end+1} = 'staticOnly'; defaults{end+1} = false;
defaults{end+1} = 'periodicPeriod'; defaults{end+1} = 2*pi;
defaults{end+1} = 'givenStaticModel'; defaults{end+1} = [];
% If zero, then the variance for the dynamical kernel for the
% rbf/matern32 etc component is not learned.
defaults{end+1} = 'learnKernelVariance'; defaults{end+1} = 0;
% Replaces optionsDyn.inverseWidth
defaults{end+1} = 'inverseWidthMult'; defaults{end+1} = 20;
defaults{end+1} = 'enableParallelism'; defaults{end+1} = false;
defaults{end+1} = 'DgtN'; defaults{end+1} = false;
% For initialising the latent space, we will account for numSharedDIms for
% the dimensionality of the shared space (that's just an initial
% calculation, the model actually decides how many shared dimensions to
% use)
defaults{end+1} = 'numSharedDims'; defaults{end+1} = 4;
defaults{end+1} = 'reconstrIters'; defaults{end+1} = 800;
% If set to true, then the actual data will be normalised immediately after
% loading them
defaults{end+1} = 'normaliseData'; defaults{end+1} = false;
defaults{end+1} = 'dataToKeep'; defaults{end+1} = -1;
% See 'utils_transformMultiLabel.m'
defaults{end+1} = 'transformMultiLabel'; defaults{end+1} = false;
defaults{end+1} = 'equalizeScales'; defaults{end+1} = false;

% This allows to rewrite the standard form of the bound F + KL as
% 2*((1-fw)*F + fw*KL), where fw is the KLweight. This means that when fw
% is 0.5 there's no difference, otherwise the KL term can be more or less
% emphasized
defaults{end+1} = 'KLweight'; defaults{end+1} = 0.5;

% adjust the "influence" of each of the partial likelihood bounds according
% to their dimensionality. Specifically, if the bound is:
% F'; defaults{end+1} = F_1 + F_2 + ... + F_M - KL, then each F_i is scaled with
% 1/dimensionality(F_i) and all F terms are rescaled back so that the
% KL part is not unbalanced (ie all the coefficients of F_i's sum to 1)
% This is an evil way of balancing the models and should be avoided!!!
% See below a better way.
defaults{end+1} = 'balanceModalities'; defaults{end+1} = false;

% Similarly to the above, this tackles the fact that sometimes
% modalities have very different number of dimensions. This field
% linearly maps the selected modalities in higher dimensions via a
% random matrix "modalityMapping" which is stored into the model. This
% matrix maps to the dimensionality of the biggest modality. If
% balanceModalityDim is a single value (true or false) this is
% replicated to every modality. If it's a cell array, it only affects
% the corresponding modalities.
% The reverse mapping can be made by simply multiplying every row of
% modality i with model.modalityMapping{i}' (this should be done e.g.
% for predictions).
defaults{end+1} = 'balanceModalityDim'; defaults{end+1} = false;
defaults{end+1} = 'optimiser'; defaults{end+1} = 'scg';
defaults{end+1} = 'saveModelDir'; defaults{end+1} = './';

% Possible options for the initX function of the models
defaults{end+1} = 'initFuncOptions'; defaults{end+1} = {};