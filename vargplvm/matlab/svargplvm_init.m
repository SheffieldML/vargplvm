% SVARGPLVM_INIT Initialise options. If the field name of 'defaults' already exists as a
% variable, the globalOpt will take this value, otherwise the default one.
%
% COPYRIGHT: Andreas C. Damianou, 2012, 2013
% VARGPLVM

if ~exist('globalOpt')
    %% Default values for the model
    defaults.experimentNo = 404;
    % Set this value to 'noSave' if you don't want to save the model
    defaults.saveName = [];
    defaults.itNo = [500 2000];
    defaults.indPoints = 100;
    defaults.latentDim = 10;
    % This should be a cell array of length M (num. of modalities) if the
    % field initial_X is set to 'custom'
    defaults.latentDimPerModel = 3;
    defaults.initVardistIters = 300;
    % How to initialise X: 
    %   'separately', means apply the initX function to each
    % of the datasets and then concatenate. 
    %   'concatenated', means first concatenate
    % the datasets and then apply the 'initX' function. 
    %   'custom' is like the "separately" options, but it implies that
    % latentDimPerModel is a cell specifying how many dimensions to use for each submodel.
    defaults.initial_X = 'separately'; % other options: 'concatenated' or 'custom'
        
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
    defaults.initX = 'ppca';
    defaults.inputScales = [];
    % Set to 1 to tie the inducing points with the latent vars. X
    defaults.fixInd = 0;
    defaults.baseKern = {'rbfardjit', 'rbfardjit'};   %{'rbfard2', 'white'};
    defaults.dynamicKern = {'rbf','white', 'bias'};
    % Like dataToKeep, but the rest indices go to the test set.
    defaults.indTr = -1;
    %%%defaults.backConstraints = 1; % Substituted with dynamicsConstrainType
    defaults.vardistCovarsMult = 2;
    defaults.dataSetNames = {};
    defaults.dataSetName = 'unknown';
    defaults.mappingInitialisation = 0;
    defaults.scale2var1 = 0;
    % Set to -1 to use all data. Set to scalar d, to only take d points from
    % each class. Set to a vector D with length equal to the number of classes C,
    % to take D(c) points from class c.
    defaults.dataPerClass = -1;
    % Set to -1 to keep all the training data, set to a number N to only keep N
    % datapoints.
    defaults.dataToKeep = -1;
    % Signal to noise ratio (initialisation for model.beta).
    defaults.initSNR = 100;
    % How many iterations to do to initialise the model with a static B. gplvm.
    % -1 means that no such initialisation will be done.
    defaults.initWithStatic = -1;
    % If initWithStatic ~= -1, this says how many iters with fixed
    % beta/sigmaf to perform.
    defaults.initWithStaticInitVardist = 300;
    % If initWithStatic ~= -1, this says what the initial SNR will be for the
    % initial static model.
    defaults.initStaticSNR = 25;
    % If true, then if also initWithStatic ~=1, we initialise the model.beta
    % and model.kern based on the initialised static model.
    defaults.initWithStaticAll = false;
    % Leave empty [] for no dynamics. 
    % Other options are: 'time', 'labels', 'data'
    defaults.dynamicsConstrainType = {};
    % If set to fals, then the dynamics kernel is not initialised with
    % bc_initdynKernel
    defaults.initDynKernel = 1;
    % A second (probably better) way to initialise the model
    defaults.initDynKernel2 = 0;
    % See bc_backConstraintsModelCreate and bc_restorePrunedModel
    defaults.labelMatrixPower = 0.5;
    % if not empty, the corresponding gradients of the kernel will be zero,
    % i.e. not learning these elements (typically for the variance of the
    % rbf/matern etc of a invcmpnd)
    defaults.fixedKernVarianceIndices = [];
    % If discrKernel is 'ones', then the discriminative kernel is
    % simply build based on a matrix with ones and minus ones, otherwise it is
    % based on a measure on the distance of each label from the mean of each
    % class.
    defaults.discrKernel = 'ones';
    % Default variance for a fixedwhite kernel
    defaults.fixedwhiteVar = 1e-5;
    % If set to some value, call it x, then after learning a constrained model,
    % (and if the function runStaticModel is called), a static model will be
    % initialised with the constrained model's latent space and learned for x iterations.
    defaults.runStaticModel = -1;
    defaults.runStaticModelInitIters = [];
    % Possible values: (none, one or both of them): 'labelsYinputs' and
    % 'labelsYoutputs'. If the first is there, then the Y of p(X|Y) is
    % augmented with the labels as C extra dimensions, where C is the total
    % number of classes (as we use 1-of-K encoding). Similarly with labelsYoutputs.
    defaults.dataConstraints = {};
    % Option to only run a static model.
    defaults.staticOnly = false;
    defaults.periodicPeriod = 2*pi;
    defaults.givenStaticModel = [];
    % If zero, then the variance for the dynamical kernel for the
    % rbf/matern32 etc component is not learned.
    defaults.learnKernelVariance = 0;
    % Replaces optionsDyn.inverseWidth
    defaults.inverseWidthMult = 20;
    defaults.enableParallelism = false;
    defaults.DgtN = false;
    % For initialising the latent space, we will account for numSharedDIms for
    % the dimensionality of the shared space (that's just an initial
    % calculation, the model actually decides how many shared dimensions to
    % use)
    defaults.numSharedDims = 4;
    defaults.reconstrIters = 800;
    % If set to true, then the actual data will be normalised immediately after
    % loading them
    defaults.normaliseData = false;
    defaults.dataToKeep = -1;
    % See 'utils_transformMultiLabel.m'
    defaults.transformMultiLabel = false;
    defaults.equalizeScales = false;
    
    % This allows to rewrite the standard form of the bound F + KL as
    % 2*((1-fw)*F + fw*KL), where fw is the KLweight. This means that when fw
    % is 0.5 there's no difference, otherwise the KL term can be more or less
    % emphasized
    defaults.KLweight = 0.5;
    
    % adjust the "influence" of each of the partial likelihood bounds according
    % to their dimensionality. Specifically, if the bound is:
    % F = F_1 + F_2 + ... + F_M - KL, then each F_i is scaled with
    % 1/dimensionality(F_i) and all F terms are rescaled back so that the
    % KL part is not unbalanced (ie all the coefficients of F_i's sum to 1)
    % This is an evil way of balancing the models and should be avoided!!!
    % See below a better way.
    defaults.balanceModalities = false;
    
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
    defaults.balanceModalityDim = false;
    defaults.optimiser = 'scg';
    defaults.saveModelDir = './';
    
    %% The final options structure will contain default values, apart for whatever
    % is already set as a variable
    fnames = fieldnames(defaults);
    for i=1:length(fnames)
        if ~exist(fnames{i}, 'var')
            globalOpt.(fnames{i}) = defaults.(fnames{i});
        else
            globalOpt.(fnames{i}) = eval(fnames{i});
        end
    end
    clear('defaults', 'fnames');
    
end