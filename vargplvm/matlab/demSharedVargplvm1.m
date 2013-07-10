% DEMSHAREDVARGPLVM1 Run the shared variational GPLVM on various kinds of data.
% DESC
%
% COPYRIGHT: Andreas C. Damianou, Carl Henrik Ek,  2011
% SEEALSO : svargplvmPrepareData
%
% VARGPLVM

%{
%---- How to run this demo:
% Vargplvm created data
close all; initLatent='ppca'; dataSetNames = {}; toyDataCreate ='vargplvm';initVardistIters = 400; itNo = [100 200 1500 500 500 400];  demSharedVargplvm1

% Fols data (make sure to use a linard2 kernel)
close all; initLatent='pca'; dataSetNames = {}; toyDataCreate ='fols';initVardistIters = 170; itNo = 300; dataToKeep=50;  demSharedVargplvm1

% Human pose data with features extracted from silhouettes
clear ;close all; initLatent='pca3'; dataSetNames = {}; toyDataCreate='humanPose';initVardistIters = 160; itNo = [200 200 200 200 200 200 200 200]; dataToKeep=418;  demSharedVargplvm1
clear ;close all; experimentNo=1;initLatent='pca3'; dataSetNames = {};toyDataCreate='humanPose';initVardistIters = 180; itNo = [500 200 200 200 200 200 200 200 200];  demSharedVargplvm1

% Human pose data with the whole silhouette (TEMPdemSharedVargplvm is like
% demSharedVargplvm but just uses more latent dimensions)
clear ;close all; experimentNo=6; imageSilhouette=1;initLatent='ppca'; latentDimPerModel=7;dataSetNames = {};toyDataCreate='humanPose';initVardistIters = 380; itNo = [500 200 200 200 200 200 200 200 200 200];  indPoints=100;TEMPdemSharedVargplvm

% TODO: ___
clear ;close all; experimentNo=5;initLatent='ppca'; mappingKern = 'rbfardjit';dataSetNames = {};toyDataCreate='humanPose';initVardistIters = 180; itNo = [500 200 200 200 200 200 200 200 200];  demSharedVargplvm1
%}

%___

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

% Define constants (in a manner that allows other scripts to parametrize
% this one).
if ~exist('experimentNo') ,  experimentNo = 404;      end
if ~exist('itNo')         ,  itNo = [500 500];              end     % Default: 2000
if ~exist('indPoints')    ,  indPoints = 80;          end     % Default: 49
if ~exist('latentDimPerModel')    ,  latentDimPerModel = 3;          end
if ~exist('latentDim'), latentDim = 5; end
if ~exist('numSharedDims'), numSharedDims = 2; end
% Set to 1 to use dynamics or to 0 to use the standard var-GPLVM
if ~exist('dynUsed')      ,  dynUsed = 0;             end
if ~exist('initVardistIters'), initVardistIters = 100;      end
if ~exist('dynamicKern')   ,     dynamicKern = {'rbf', 'white', 'bias'}; end
if ~exist('mappingKern')   ,  mappingKern = {'rbfard2', 'white'}; end
% if ~exist('mappingKern'),  mappingKern = {'rbfard2', 'bias', 'white'}; end
% Set to 1 to tie the inducing points with the latent vars. X
if ~exist('fixInd')        ,     fixInd = 0;    end
% 0.1 gives around 0.5 init.covars. 1.3 biases towards 0.
if ~exist('vardistCovarsMult'),  vardistCovarsMult=1.3;                  end
% Set to empty value {} to work with toy data
if ~exist('dataSetNames')    ,    dataSetNames = {};    end
if ~exist('invWidthMultDyn'),    invWidthMultDyn = 100;                     end
if ~exist('invWidthMult'),       invWidthMult = 5;                     end
if ~exist('initX'),     initX ='ppca';   end % That's for the dynamics initialisation
if ~exist('initLatent'),     initLatent ='ppca';   end % That's for the whole latent space init.
if ~exist('dataToKeep'), dataToKeep = -1; end
if ~exist('toyDataCreate'), toyDataCreate = 'fols'; end
if ~exist('doPredictions'), doPredictions = 0; end
if ~exist('dataType'), dataType = 'default'; end
if ~exist('enableParallelism'), enableParallelism = 1; end

if ~exist('Yall')
    svargplvmPrepareData;
end

%%

%-- Load datasets
for i=1:numberOfDatasets
    Y = Yall{i};
    if dataToKeep ~= -1 && dataToKeep <= size(Y,1)
        Y = Y(1:dataToKeep,:);
    end
    dims{i} = size(Y,2);
    N{i} = size(Y,1);
    indTr = 1:N{i};
    indTs = setdiff(size(Y,1), indTr);
    Ytr{i} = Y(indTr,:); %Yts = Y(indTs,:);
    t{i} = linspace(0, 2*pi, size(Y, 1)+1)'; t{i} = t{i}(1:end-1, 1);
    timeStampsTraining{i} = t{i}(indTr,1); %timeStampsTest = t(indTs,1);
    d{i} = size(Ytr{i}, 2);
end

for i=2:numberOfDatasets
    if N{i} ~= N{i-1}
        error('The number of observations in each dataset must be the same!');
    end
end

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
end

X_init = svargplvmInitLatentSpace(Ytr, d, options, initLatent, latentDim, latentDimPerModel, numSharedDims);

latentDim = size(X_init,2);

% Free up some memory
clear('Y')



%-- Create the sub-models: Assume that for each dataset we have one model.
% This can be changed later, as long as we find a reasonable way to
% initialise the latent spaces.
for i=1:numberOfDatasets
    %---- Here put some code to assign X to the global common X which must
    % be created by doing pca in the concatenation of Y's...After this
    % point, model{i}.X will be the same for all i's. TODO...
    fprintf(1,'# Creating the model...\n');
    options{i}.initX = X_init;
    model{i} = vargplvmCreate(latentDim, d{i}, Ytr{i}, options{i});
    model{i}.X = X_init; %%%%%%%
    model{i} = vargplvmParamInit(model{i}, model{i}.m, model{i}.X);
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
    
    %-------- Add dynamics to the model -----
    if dynUsed
        fprintf(1,'# Adding dynamics to the model...\n');
        optionsDyn{i}.type = 'vargpTime';
        optionsDyn{i}.t=timeStampsTraining{i};
        optionsDyn{i}.inverseWidth=invWidthMultDyn; % Default: 100
        optionsDyn{i}.initX = X_init; % initX; % That's for the dynamics
        
        kern = kernCreate(t{i}, dynamicKern); % Default: {'rbf','white','bias'}
        
        %-- Initialize each element of the compound kernel
        svargplvmInitDynKern
        
        optionsDyn{i}.kern = kern;
        optionsDyn{i}.vardistCovars = vardistCovarsMult; % 0.23 gives true vardist.covars around 0.5 (DEFAULT: 0.23) for the ocean dataset
        
        % Fill in with default values whatever is not already set
        optionsDyn{i} = vargplvmOptionsDyn(optionsDyn{i});
        model{i} = vargplvmAddDynamics(model{i}, 'vargpTime', optionsDyn{i}, optionsDyn{i}.t, 0, 0,optionsDyn{i}.seq);
        
        fprintf(1,'# Further calibration of the initial values...\n');
        model{i} = vargplvmInitDynamics(model{i},optionsDyn{i});
        
        %  to also not learn the last kernel's variance
        if numel(kern.comp) > 1 && exist('learnSecondVariance') && ~learnSecondVariance
            fprintf(1,'# The variance for %s in the dynamics is not learned!\n',kern.comp{end}.type)
            model{i}.dynamics.learnSecondVariance = 0;
            model{i}.dynamics.kern.comp{end}.inverseWidth = model{i}.dynamics.kern.comp{1}.inverseWidth/10; %%% TEMP
        end
    end
    
    model{i}.beta=1/(0.01*var(model{i}.m(:)));
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
capName = dataType;
capName(1) = upper(capName(1));
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

%%


display = 1;
%%%% Optimisation
% do not learn beta and sigma_f for few iterations for intitialization
if initVardistIters ~=0
    model.initVardist = 1;
    model = svargplvmPropagateField(model,'initVardist', model.initVardist);
    fprintf(1,'# Intitiliazing the variational distribution...\n');
    model = svargplvmOptimise(model, display, initVardistIters); % Default: 20
    %fprintf(1,'1/b = %.4d\n',1/model.beta);
    
    modelInitVardist = model;
    model.initVardistIters=initVardistIters;
end

model.initVardist = 0;
model = svargplvmPropagateField(model,'initVardist', model.initVardist);

% Optimise the model.
model.iters = 0;
for i=1:length(itNo)
    iters = itNo(i); % default: 2000
    fprintf(1,'\n# Optimising the model for %d iterations (session %d)...\n',iters,i);
    model = svargplvmOptimise(model, display, iters);
    model.iters = model.iters + iters;
    % Save model
    fprintf(1,'# Saving %s\n',fileToSave);
    save(fileToSave, 'model', 'prunedModelInit');
end

for i=1:numberOfDatasets
    figure, bar(model.comp{i}.kern.comp{1}.inputScales);
    title(['Final scales for dataset ' num2str(i)]);
end

if toyData && strcmp(toyDataCreate,'fols')
    figure,plot(model.X(:,1)), hold on, plot(model.X(:,2),'r')
    hold on, plot(model.X(:,3),'g')
end

%%
%--------------------- PREDICTIONS --------------%
if ~doPredictions
    return
end
if model.numModels ~=2
    error('Predictions cannot be made (currently) with sharedVargplvm unless the number of submodels used is 2.');
end
%%

obsMod = 1; % one of the involved sub-models (the one for which we have the data)
infMod = setdiff(1:2, obsMod);

% Find the dimensions that are shared for obsMod and infMod
if ~exist('sharedDims')
    thresh = max(model.comp{obsMod}.kern.comp{1}.inputScales) * 0.001;
    retainedScales{obsMod} = find(model.comp{obsMod}.kern.comp{1}.inputScales > thresh);
    thresh = max(model.comp{infMod}.kern.comp{1}.inputScales) * 0.001;
    retainedScales{infMod} = find(model.comp{infMod}.kern.comp{1}.inputScales > thresh);
    sharedDims = intersect(retainedScales{obsMod}, retainedScales{infMod});
end

% Find X_* only for the shared dimensions (Xs*):
if ~exist('privateDims')
    privateDims{infMod} = setdiff(1:model.comp{obsMod}.q, sharedDims);
end

testOnTraining=0;

if testOnTraining
    i = size(model.comp{obsMod}.y,1); % last observation
    % Find p(X_* | Y_*): If Y* is not taken from the tr. data, then this step
    % must be an optimisation step of q(x*). X_* and Y_* here refer to the
    % spaces of the submodel obsMod.
    y_star = model.comp{obsMod}.y(i);
    x_star = model.comp{obsMod}.vardist.means(i,:);
    varx_star = model.comp{obsMod}.vardist.covars(i,:);
    [ind, distInd] = nn_class(model.X(:,sharedDims), x_star(:,sharedDims), 10, 'euclidean');
    
    ZpredMu = zeros(length(ind), size(Ytoy{infMod},2));
    ZpredSigma = zeros(length(ind), size(Ytoy{infMod},2));
    for k=1:length(ind)
        [ZpredMu(k,:), ZpredSigma(k,:)] = vargplvmPosteriorMeanVar(model.comp{infMod}, model.X(ind(k),:));
    end
else
    testInd = 25;
    Yts = Y_test(testInd,:); % Now this is a matrix
    for i=1:size(Yts,1)
        % initialize the latent points using the nearest neighbour
        % from the training data
        dst = dist2(Yts(i,:), model.comp{obsMod}.y);
        [mind, mini] = min(dst);
        
        Init(i,:) = model.vardist.means(mini,:);
        vardistx = vardistCreate(model.comp{obsMod}.vardist.means(mini,:), model.q, 'gaussian');
        vardistx.covars = model.comp{obsMod}.vardist.covars(mini,:);
        model.comp{obsMod}.vardistx = vardistx;
        display=1;
        iters = 100;
        % Find p(X_* | Y_*) which is approximated by q(X_*)
        [x_star, varx_star, modelUpdated] = vargplvmOptimisePoint(model.comp{obsMod}, vardistx, Yts(i,:), display, iters);%%%
        numberOfNN = 10;
        % Now we selected a datapoint X_* by taking into account only the
        % private dimensions for Y. Now, based on the shared dimensions of
        % that, we select the closest (in a NN manner) X from the training data.
        [ind, distInd] = nn_class(model.X(:,sharedDims), x_star(:,sharedDims), numberOfNN, 'euclidean');
        
        x = model.X(ind,:);
        % Careful:  Since x is coming from the training set, it is somewhat
        % natural that the predicted y's will also closely match the ones
        % from the training set. What makes it interesting, is to show that
        % the shared dimensions are the ones controlling the common
        % ambiguities, so, actually it really matters what indices of
        % training X's we find with the NN based on sharedDims.
        mu = vargplvmPosteriorMeanVar(model.comp{infMod}, x);
        
        ZpredK{i} = mu;
        
        ZpredMu(i,:) = mu(1);
        
        %ZpredSigma(i,:) = sigma;
    end
    
    errsumFull = sum((ZpredMu - Z_test(testInd,:)).^2);
    errorFull = mean(errsumFull);
    
    if toyData && strcmp(toyDataCreate,'humanPose')
        for j = 1:1:length(testInd)
            handle = xyzankurVisualise(ZpredK{j}(1,:),1);
            for k = 2:length(ind)
                %  xyzankurModify(handle,ZpredK{j}(k,:));            handle = xyzankurVisualise(ZpredK{j}(1,:),1);
                handle = xyzankurVisualise(ZpredK{j}(k,:),k);
                title(['Pose'])
                pause;
            end
        end
    end
end

