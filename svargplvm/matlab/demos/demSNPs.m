% DEMFMRISVARGPLVM2 Run the shared variational GPLVM on fmri data
%
% COPYRIGHT: Andreas C. Damianou
%
% SHEFFIELDML


%
% clear;experimentNo = 1; mappingKern = {'linard2','white','bias'}; indPoints = 120;initVardistITers = 200; latentDimPerModel=20;demSNPs % Good!
% clear;experimentNo = 2; initVardistIters = 220; indPoints = 60; demSNPs %(lovelace)
% clear;experimentNo = 3; mappingKern = {'linard2','white','bias'}; indPoints = 120;  latentDimPerModel=20;initial_X='together';demSNPs %(laplace)
% clear;experimentNo = 4; initVardistIters = 200; indPoints = 60; initial_X='together';demSNPs %(laplace)

% clear;experimentNo = 7; mappingKern = {'linard2','white','bias'}; initial_X = 'separately'; indPoints = 109; latentDimPerModel=20; initVardistIters = 200; demSNPs;

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

% Define constants (in a manner that allows other scripts to parametrize
% this one).
if ~exist('experimentNo') ,  experimentNo = 404 ;      end
if ~exist('itNo')         ,  itNo = [1000 500 500 500];              end     % Default: 2000
if ~exist('indPoints')    ,  indPoints = 60;          end     % Default: 49
if ~exist('latentDimPerModel')    ,  latentDimPerModel = 12;          end
% Set to 1 to use dynamics or to 0 to use the standard var-GPLVM
if ~exist('initVardistIters'), initVardistIters = 120;      end
%if ~exist('mappingKern')   ,  mappingKern = {'rbfard2', 'white'}; end
if ~exist('mappingKern')   ,  mappingKern = 'rbfardjit'; end
% Set to 1 to tie the inducing points with the latent vars. X
if ~exist('fixInd')        ,     fixInd = 0;    end
% Set to empty value {} to work with toy data
if ~exist('dataSetNames')    ,    dataSetNames = {};    end
if ~exist('initLatent'),     initLatent ='ppca';   end % That's for the whole latent space init.
if ~exist('dataType'), dataType = 'snps'; end
if ~exist('invWidthMult'),       invWidthMult = 5;                     end
if ~exist('enableParallelism'), enableParallelism = 0; end
if ~exist('enableDgtN'), enableDgtN = true; end
if ~exist('initial_X'), initial_X = 'separately'; end
if ~exist('initSNR'), initSNR = 100; end
if ~exist('associationMode'), associationMode = 'eth-glu'; end
%%



%%
Y = svargplvmLoadData('bio/smith_yeast');

switch associationMode
    case 'snps-expr-ethanol'
        Ytr{1} = Y.expr_ethanol;
        Ytr{2} = Y.snps_ethanol;
        dataSetNames = {'expr_ethanol','snps_ethanol'};
    case 'snps-expr-glucose'
        Ytr{1} = Y.expr_glucose;
        Ytr{2} = Y.snps_glucose;
        dataSetNames = {'expr_glucose','snps_glucose'};
    case 'snps-expr-combined'
        Ytr{1} = [Y.expr_ethanol Y.expr_glucose];
        Ytr{2} = [Y.snps_ethanol Y.snps_glucose];
        dataSetNames = {'expr_ethanol_glucose','snps_ethanol_glucoe'};
    otherwise % 'eth-glu'
        Ytr{1} = Y.expr_ethanol;
        Ytr{2} = Y.expr_glucose;
        dataSetNames = {'expr_ethanol','expr_glucose'};
end
clear('Y')

%{
%--- center and scale data
for i=1:length(Ytr)
    bias = mean(Ytr{i});
    scale = std(Ytr{i});
    scale(find(scale==0)) = 1;
    % Remove bias and apply scale.
    for j = 1:size(Ytr{i},2)
        Ytr{i}(:, j) = Ytr{i}(:, j) - bias(j);
        if scale(j)
            Ytr{i}(:, j) = Ytr{i}(:, j)/scale(j);
        end
    end
end
%--
%}


numberOfDatasets = length(Ytr);

%-- Load datasets
for i=1:numberOfDatasets
    Y = Ytr{i};
    dims{i} = size(Y,2);
    N{i} = size(Y,1);
    d{i} = size(Ytr{i}, 2);
end

for i=2:numberOfDatasets
    if N{i} ~= N{i-1}
        error('The number of observations in each dataset must be the same!');
    end
end

%% 
%-- Options for the models
for i=1:numberOfDatasets
    % Set up models
    options{i} = vargplvmOptions('dtcvar');
    options{i}.kern = mappingKern; %{'rbfard2', 'bias', 'white'};
    %indPoints = 80; %%%%%
    options{i}.numActive = indPoints;
    options{i}.optimiser = 'scg2';
    
    % !!!!! Be careful to use the same type of scaling and bias for all
    % models!!!
    
    % scale = std(Ytr);
    % scale(find(scale==0)) = 1;
    %options.scaleVal = mean(std(Ytr));
    options{i}.scaleVal = sqrt(var(Ytr{i}(:)));
    if fixInd
        options{i}.fixInducing=1;
        options{i}.fixIndices=1:size(Ytr{i},1);
    end
    options{i}.enableDgtN = enableDgtN;
end



mAll=[];
%-- Create the normalised version of the datasets and concatenate
%!!!! Doesn't work if Y{i}'s have different sizes!!
for i=1:length(Ytr)
    % Compute m, the normalised version of Ytr (to be used for
    % initialisation of X)
    bias = mean(Ytr{i});
    scale = ones(1, d{i});
    
    if(isfield(options{i},'scale2var1'))
        if(options{i}.scale2var1)
            scale = std(Ytr{i});
            scale(find(scale==0)) = 1;
            if(isfield(options{i}, 'scaleVal'))
                warning('Both scale2var1 and scaleVal set for GP');
            end
        end
    end
    if(isfield(options{i}, 'scaleVal'))
        scale = repmat(options{i}.scaleVal, 1, d{i});
    end
    
    % Remove bias and apply scale.
    m{i} = Ytr{i};
    for j = 1:d{i}
        m{i}(:, j) = m{i}(:, j) - bias(j);
        if scale(j)
            m{i}(:, j) = m{i}(:, j)/scale(j);
        end
    end
end
% Clear some variables
clear('Y','bias','scale','ind2');

if strcmp(initial_X, 'separately')
    fprintf('# Initialising X by performing ppca in each observed (scaled) dataset separately and then concatenating...\n');
    X_init{1} = ppcaEmbed(m{1},latentDimPerModel);
    X_init{2} = ppcaEmbed(m{2},latentDimPerModel);
    X_init = [X_init{1} X_init{2}];
else
    fprintf('# Initialising X by performing ppca in concatenated observed (scaled) data...\n');
    X_init = ppcaEmbed([m{1} m{2}], latentDimPerModel*2);
end
%-----------------

latentDim = size(X_init,2);





%-- Create the sub-models: Assume that for each dataset we have one model.
% This can be changed later, as long as we find a reasonable way to
% initialise the latent spaces.
for i=1:numberOfDatasets
    %---- Here put some code to assign X to the global common X which must
    % be created by doing pca in the concatenation of Y's...After this
    % point, model{i}.X will be the same for all i's. TODO...
    fprintf(1,'# Creating the model...\n');
    options{i}.initX = X_init;
    %%%%%%!!!!!!! TEMP: There is a function which overrides the default...
    %model{i} = TEMPvargplvmCreate(latentDim, d{i}, Ytr{i}, options{i},0);
    model{i} = vargplvmCreate(latentDim, d{i}, Ytr{i}, options{i});
    model{i}.X = X_init; %%%%%%%
    model{i} = vargplvmParamInit(model{i}, m{i}, model{i}.X);
    model{i}.X = X_init; %%%%%%%
    
    inpScales = invWidthMult./(((max(model{i}.X)-min(model{i}.X))).^2); % Default 5
    %inpScales(:) = max(inpScales); % Optional!!!!!
    model{i}.kern.comp{1}.inputScales = inpScales;
    
    if strcmp(model{i}.kern.type, 'rbfardjit')
        model{i}.kern.inputScales = model{i}.kern.comp{1}.inputScales;
    end
    params = vargplvmExtractParam(model{i});
    model{i} = vargplvmExpandParam(model{i}, params);
    model{i}.vardist.covars = 0.5*ones(size(model{i}.vardist.covars)) + 0.001*randn(size(model{i}.vardist.covars));
    
    
    if model{i}.DgtN
        model{i}.beta=initSNR/var(model{i}.mOrig(:));
    else
        model{i}.beta=initSNR/var(model{i}.m(:));
    end
    prunedModelInit{i} = vargplvmPruneModel(model{i});
    %disp(model{i}.vardist.covars)
end

if experimentNo == -1
    experimentNo = globalExperimentNumber('sharedVargplvm', 'dataType');
end

%modelInit = model;%%%TEMP

%--  Unify models into a structure
model = svargplvmModelCreate(model);
model.dataSetNames = dataSetNames;
model.initLatent = initLatent;
model.experimentNo = experimentNo;
model.dataType = dataType;
%%---
if ~isempty(dataType)
    capName = dataType;
    capName(1) = upper(capName(1));
else
    capName = [];
end
modelType = model.type;
modelType(1) = upper(modelType(1));
fileToSave = ['dem' capName modelType num2str(experimentNo) '.mat'];
%%---


%-- Define what level of parallelism to use (w.r.t submodels or/and w.r.t
% datapoints).
if enableParallelism
    if model.numModels > 8
        fprintf('# Parallel computations w.r.t the submodels!\n');
        model.parallel = 1;
        model = svargplvmPropagateField(model,'parallel', 1);
    elseif model.N > 15
        fprintf('# Parallel computations w.r.t the datapoints!\n');
        model.vardist.parallel = 1;
        for i=1:model.numModels
            model.comp{i}.vardist.parallel = 1;
        end
    end
end

% Force kernel computations
params = svargplvmExtractParam(model);
model = svargplvmExpandParam(model, params);
model.initial_X = initial_X;
%%


display = 1;
%%%% Optimisation
% do not learn beta and sigma_f for few iterations for intitialization
for i=1:length(initVardistIters)
    if ~initVardistIters(i) > 0
        continue
    end
    model.initVardist = 1; model.learnSigmaf = 0;
    model = svargplvmPropagateField(model,'initVardist', model.initVardist);
    model = svargplvmPropagateField(model,'learnSigmaf', model.learnSigmaf);
    fprintf(1,'# Intitiliazing the variational distribution for %d iters...\n', initVardistIters(i));
    model = svargplvmOptimise(model, display, initVardistIters(i)); % Default: 20
    %fprintf(1,'1/b = %.4d\n',1/model.beta);
    
    model.initVardistIters=initVardistIters(i);
    
    prunedModel = svargplvmPruneModel(model);
    fprintf(1,'# Saving %s\n',fileToSave);
    save(fileToSave, 'prunedModel', 'prunedModelInit');
end
modelInitVardist = model;

    
model.initVardist = 0; model.learnSigmaf=1;
model = svargplvmPropagateField(model,'initVardist', model.initVardist);
model = svargplvmPropagateField(model,'learnSigmaf', model.learnSigmaf);

% Optimise the model.
model.iters = 0;
for i=1:length(itNo)
    iters = itNo(i); % default: 2000
    fprintf(1,'\n# Optimising the model for %d iterations (session %d)...\n',iters,i);
    model = svargplvmOptimise(model, display, iters);
    model.iters = model.iters + iters;
    % Save model
    prunedModel = svargplvmPruneModel(model);
    fprintf(1,'# Saving %s\n',fileToSave);
    save(fileToSave, 'prunedModel', 'prunedModelInit');
end

svargplvmShowScales(model,false)
SNR = svargplvmSNR(model)

