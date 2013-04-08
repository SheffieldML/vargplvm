% DEMFMRISVARGPLVM3 Run the shared variational GPLVM on fmri data
% 
% Similar to demFmriSvargplvm2 but here we also add the labels as an extra
% modality
%
% COPYRIGHT: Andreas C. Damianou
%
% SHEFFIELDML


% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

% Define constants (in a manner that allows other scripts to parametrize
% this one).
if ~exist('experimentNo') ,  experimentNo = 200 ;      end
if ~exist('itNo')         ,  itNo = [700 1500 500 500];              end     % Default: 2000
if ~exist('indPoints')    ,  indPoints = 80;          end     % Default: 49
if ~exist('latentDims')    ,  latentDims = 25;          end
% Set to 1 to use dynamics or to 0 to use the standard var-GPLVM
if ~exist('dynUsed')      ,  dynUsed = 0;             end
if ~exist('initVardistIters'), initVardistIters = 500;      end
if ~exist('dynamicKern')   ,     dynamicKern = {'rbf', 'white', 'bias'}; end
%if ~exist('mappingKern')   ,  mappingKern = {'rbfard2', 'white'}; end
if ~exist('mappingKern')   ,  mappingKern = 'rbfardjit'; end
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
% Initialise latent space by applying @initLatent separately and then
% concatenate, or apply @initX on the already concatenated output space?
% Options: 'together', 'separately'
if ~exist('initial_X'), initial_X = 'together'; end 
if ~exist('dataType'), dataType = 'fmriWithLbls'; end
if ~exist('enableParallelism'), enableParallelism = 1; end
if ~exist('doPredictions'), doPredictions = 0; end
if ~exist('keepDiary'), keepDiary = false; end
%%
disp('# Loading and preparing datasets...')

totalRuns900 = 6;
totalRuns400 = 7;

%-- Training data, 1st dataset
%load fmri900Processed2
load fmri900Processed4
dataToKeep = find(runNumber <= totalRuns900);
Yall{1} = Y(dataToKeep,:);
activityAll{1} = activity(dataToKeep);
runNumberAll{1} = runNumber(dataToKeep);

% Test data, 1st dataset
dataToKeep = find(runNumber > totalRuns900);
YallTest{1} = Y(dataToKeep,:);
activityTest{1} = activity(dataToKeep);
runNumberTest{1} = runNumber(dataToKeep);

%-- Training data, 2nd dataset
%load fmri400Processed2
load fmri400Processed4
dataToKeep = find(runNumber <= totalRuns400);
Yall{2} = Y(dataToKeep,:);
activityAll{2} = activity(dataToKeep);
runNumberAll{2} = runNumber(dataToKeep);

% Test data, 2nd dataset
dataToKeep = find(runNumber > totalRuns400);
YallTest{2} = Y(dataToKeep,:);
activityTest{2} = activity(dataToKeep);
runNumberTest{2} = runNumber(dataToKeep);

%--- Make training datasets compatible: Each of the Y{i}, activity{i}, runNumber{i}
% will correspond to the same object. This is done by rearranging these
% original cell-arrays. For some elements, one of the two datasets has one
% more entry. We remove this entry from the corresponding dataset.

clear('Y','runNumber','activity')

j=1;
allElements = unique(activityAll{1});
for i=1:length(allElements)
    ind1 = find(strcmp(allElements{i},activityAll{1}));
    ind2 = find(strcmp(allElements{i},activityAll{2}));
    l1 = length(ind1);
    l2 = length(ind2);
    if l2 < l1
        ind1 = ind1(1:l2);
    elseif l1 < l2
        ind2 = ind2(1:l1);
    end
    l = min(l1,l2);
    Y{1}(j:j+l-1,:) = Yall{1}(ind1,:);
    activity{1}(j:j+l-1) = activityAll{1}(ind1);
    runNumber{1}(j:j+l-1) = runNumberAll{1}(ind1);
    Y{2}(j:j+l-1,:) = Yall{2}(ind2,:);
    activity{2}(j:j+l-1) = activityAll{2}(ind2);
    runNumber{2}(j:j+l-1) = runNumberAll{2}(ind2)';
    j = j+l;
end
activity{1} = activity{1}'; activity{2} = activity{2}';
runNumber{1} = runNumber{1}'; runNumber{2} = runNumber{2}';
Yall = Y; clear('Y','activityAll','runNumberAll');
numberOfDatasets = 2;


% Remove entries from test sets that are nonexistent in the training set
% (e.g. objects or animals that don't exist in tr. set).
% TODO...

%%

%-- Load datasets
perm = randperm(size(Yall{1},1));
Ntr = size(Yall{1},1)./1.5;
activityAll = activity; clear('activity');
if strcmp([activityAll{1}{:}], [activityAll{2}{:}]) ~= 1
    error('Activity matrices must be the same')
end

for i=1:numberOfDatasets
    Y = Yall{i};
    dims{i} = size(Y,2);
    N{i} = Ntr;
    if N{i} ~= size(Y,1)
        indTr = perm(1:N{i});
    else
        indTr = 1:N{i};
    end
    indTs = setdiff(size(Y,1), indTr);
    indTr = sort(indTr); indTs = sort(indTs);
    Ytr{i} = Y(indTr,:); %Yts = Y(indTs,:);
    t{i} = linspace(0, 2*pi, size(Y, 1)+1)'; t{i} = t{i}(1:end-1, 1);
    timeStampsTraining{i} = t{i}(indTr,1); %timeStampsTest = t(indTs,1);
    d{i} = size(Ytr{i}, 2);
    if i==1
        activity = activityAll{i}(indTr,:);
        activityTs = activityAll{i}(indTs,:);
    end
    Yts{i} = Y(indTs,:);
end



clear('activityAll');

restIndTr = find(strcmp(activity,'rest'));
animalIndTr = find(strncmp(activity,'a-',2));
otherIndTr = find(strncmp(activity,'t-',2) | strncmp(activity,'o-',2));
if length(unique([restIndTr; otherIndTr; animalIndTr])) ~= length(activity)
    error('Something is wrong with the indices')
end

for i=2:numberOfDatasets
    if N{i} ~= N{i-1}
        error('The number of observations in each dataset must be the same!');
    end
end
N = N{1};

numClasses = 3;
Ytr{end+1} = zeros(N,numClasses);
% 1-of K encoding: animal is 100, other is 010, rest is 001
Ytr{end}(animalIndTr,1) = 1;
Ytr{end}(otherIndTr,2) = 1;
Ytr{end}(restIndTr,3) = 1;
d{end+1} = numClasses;
numberOfDatasets = numberOfDatasets+1;

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



initFunc = str2func([initLatent 'Embed']);
switch initial_X
    case 'separately'
        fprintf('# Initialising the latent space by applying %s to each output space and then concatenating...\n', initLatent)
        X_init{1} = initFunc(m{1},round(latentDims/2));
        X_init{2} = initFunc(m{2},round(latentDims/2));
        X_init{3} = m{3}; % labels
        X_init = [X_init{1} X_init{2} X_init{3}];
    case 'together'
        fprintf('# Initialising the latent space by applying %s to the concatenated output spaces...\n', initLatent)
        X_init = initFunc([m{1} m{2} m{3}], latentDims);
end

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
    options{i}.initSNR = 100;
    %%%%%%!!!!!!! TEMP: There is a function which overrides the default...
    %model{i} = TEMPvargplvmCreate(latentDim, d{i}, Ytr{i}, options{i},0);
    model{i} = vargplvmCreate(latentDim, d{i}, Ytr{i}, options{i});
    model{i}.X = X_init; %%%%%%%
    if isfield(model{i}, 'mOrig')
        model{i} = vargplvmParamInit(model{i}, model{i}.mOrig, model{i}.X, options{i});
    else
        model{i} = vargplvmParamInit(model{i}, model{i}.m, model{i}.X, options{i});
    end
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
    
    %model{i}.beta=1/(0.01*var(model{i}.y(:))); %paramInit takes care of this
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

if keepDiary
    diary(['LOG_' fileToSave(1:end-3) '.txt']);
end
%%

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
    model.initVardist = 1; model.learnSigmaf = 0;
    model = svargplvmPropagateField(model,'initVardist', model.initVardist);
    model = svargplvmPropagateField(model,'learnSigmaf', model.learnSigmaf);
    fprintf(1,'# Intitiliazing the variational distribution for %d iters...\n', initVardistIters);
    model = svargplvmOptimise(model, display, initVardistIters); % Default: 20
    %fprintf(1,'1/b = %.4d\n',1/model.beta);
    
    modelInitVardist = model;
    model.initVardistIters=initVardistIters;
end
prunedModelInit = svargplvmPruneModel(model);
fprintf(1,'# Saving %s\n',fileToSave);
save(fileToSave, 'prunedModelInit', 'activity','runNumber');

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
    save(fileToSave, 'prunedModel', 'prunedModelInit', 'activity','runNumber');
end

for i=1:numberOfDatasets
    figure, bar(model.comp{i}.kern.comp{1}.inputScales);
    title(['Final scales for dataset ' num2str(i)]);
end



%%
%---------------------------- PREDICTIONS ---------------
if ~doPredictions
    return
end

mNew{1} = model.comp{1}.m;
mNew{2} = model.comp{2}.m;
model.comp{1}.m = model.comp{1}.mOrig;
model.comp{2}.m = model.comp{2}.mOrig;

obsMod = 1; % one of the involved sub-models (the one for which we have the data)
infMod = setdiff(1:2, obsMod);

% Find the dimensions that are shared for obsMod and infMod
if ~exist('sharedDims')
    s1 = model.comp{obsMod}.kern.comp{1}.inputScales;
    s2 = model.comp{infMod}.kern.comp{1}.inputScales;
    % Normalise values between 0 and 1
    s1 = s1 / max(s1);
    s2 = s2 / max(s2);
    
    %  thresh = max(model.comp{obsMod}.kern.comp{1}.inputScales) * 0.001;
    thresh = 0.005;
    
    retainedScales{obsMod} = find(s1 > thresh);
    %thresh = max(model.comp{infMod}.kern.comp{1}.inputScales) * 0.001;
    retainedScales{infMod} = find(s2  > thresh);
    sharedDims = intersect(retainedScales{obsMod}, retainedScales{infMod});
end

% Find X_* only for the shared dimensions (Xs*):
if ~exist('privateDims')
    privateDims = setdiff(1:model.comp{obsMod}.q, sharedDims);
end

Nstar = 30;
if ~exist('testOnTraining'),  testOnTraining=1; end

if testOnTraining
    testInd = randperm(size(Ytr{1},1));
    testInd = testInd(1:Nstar);
else
    testInd = randperm(size(YallTest{obsMod},1));
    testInd = testInd(1:Nstar);
    Yts = YallTest{obsMod}(testInd,:); 
end

%%


currentLabel = [];
for i=1:length(testInd)
    curInd = testInd(i);
    if testOnTraining
        currentLabel{i} = activity(curInd);
        fprintf('\n# Original: %s\n', cell2mat(currentLabel{i}));
        % Find p(X_* | Y_*): If Y_* comes from the training data, then this
        % is simply the corresponding indices of q(X).
        y_star = model.comp{obsMod}.y(curInd);
        x_star = model.comp{obsMod}.vardist.means(curInd,:);
        varx_star = model.comp{obsMod}.vardist.covars(curInd,:);
    else
        currentLabel{i} = activityTest{obsMod}(testInd(curInd));
        fprintf('\n# Original: %s\n', cell2mat(currentLabel{i}));
        % Find p(X_* | Y_*): If Y* is not taken from the tr. data, then this step
        % must be an optimisation step of q(x*). X_* and Y_* here refer to the
        % spaces of the submodel obsMod.
        % Initialize the latent point using a point of the training data which has the same label
        % OPTION 1 ---
        mini = find(strcmp(activity{obsMod},currentLabel{i}));
        mini = mini(1);
        % OPTION 2 ---
        % initialize the latent points using the nearest neighbour
        % from the training data
        % dst = dist2(Yts(i,:), model.comp{obsMod}.y);
        % [mind, mini] = min(dst);
        %---
        Init(i,:) = model.vardist.means(mini,:);
        fprintf('# Init. vardist with y with label: %s\n', cell2mat(activity{obsMod}(mini)));
        vardistx = vardistCreate(model.comp{obsMod}.vardist.means(mini,:), model.q, 'gaussian');
        vardistx.covars = model.comp{obsMod}.vardist.covars(mini,:);
        model.comp{obsMod}.vardistx = vardistx;
        display=1;
        iters = 3;
        % Find p(X_* | Y_*) which is approximated by q(X_*)
        [x_star, varx_star, modelUpdated] = vargplvmOptimisePoint(model.comp{obsMod}, vardistx, Yts(i,:), display, iters);%%%    
    end
    % Now we found a datapoint X_* by taking into account all
    % dimensions. Now, based only on the shared dimensions of
    % that, we select the closest (in a NN manner) X from the training
    % data. If we test on training data,the first X found here is of course
    % going to be the given x, so this is unsurprising. The interesting
    % thing is to see the 2nd,3rd...ranked X's.
    numberNN = 5;
    % !!! Correct for inference: This tests if the shared dimensions controls
    % indeed the classes
     [ind, distInd] = nn_class(model.X(:,sharedDims), x_star(:,sharedDims),numberNN, 'euclidean');
    % %%%!!!TEMP: this tests if the model can at least classify X's
    % [ind, distInd] = nn_class(model.X(:,retainedScales{obsMod}), x_star(:,retainedScales{obsMod}), numberNN, 'euclidean');
    
    ZpredMu = zeros(length(ind), size(model.comp{infMod}.y,2));
    ZpredSigma = zeros(length(ind), size(model.comp{infMod}.y,2));
    indicesX{i} = [];
    fprintf('# NN labels for X:\n');
    % Find p(y_*|x_*) for every x_*
    % Careful:  Since x is coming from the training set, it is somewhat
    % natural that the predicted y's will also closely match the ones
    % from the training set. What makes it interesting, is to show that
    % the shared dimensions are the ones controlling the common
    % ambiguities, so, actually it really matters what indices of
    % training X's we find with the NN based on sharedDims.
    for k=1:length(ind)
        [ZpredMu(k,:), ZpredSigma(k,:)] = vargplvmPosteriorMeanVar(model.comp{infMod}, model.X(ind(k),:));
        indicesX{i} = [indicesX{i} ',' activity(ind(k))];
        fprintf('   %s\n',cell2mat(activity(ind(k))));
    end
    indicesX{i} = indicesX{i}(2:end);
    
    indYnn = [];
    % Do a NN on the DATA space, for every predicted output.
    for j=1:size(ZpredMu,1)
        [indYnn(j), distInd] = nn_class(model.comp{infMod}.y, ZpredMu(j,:),1,'euclidean');
    end
    
    fprintf('# Predicted:\n');
    Results{i} = [];
    for j=1:length(indYnn)
        Results{i} = [Results{i} ',' activity(indYnn(j))];
        fprintf('  %s\n',cell2mat(activity(indYnn(j))));
    end
    Results{i} = Results{i}(2:end);
end

%% ---- Now summarize
classes = {'a','r','o'};

% Use the following to not use a sample prior on classes
priorOnClasses = [sum(strncmp(activity,'a-',2)), sum(strcmp(activity,'rest')), sum(strncmp(activity,'o-',2) | strncmp(activity,'t-',2))];
for i=1:length(priorOnClasses)
    priorOnClasses(i) = priorOnClasses(i)/length(activity);
end

% Use the following to not use a uniform prior on classes
%priorOnClasses = repmat(1/length(classes), 1, length(classes));


numTop = 1; weights = 1; 
%numTop = 2; weights = [.66, .33];

for i=1:length(classes)
    accuracyMRD.(classes{i}) = 0; accuracyMRD.([classes{i} '_N']) = zeros(1, numTop);
    accuracyRAND.(classes{i}) = 0; accuracyRAND.([classes{i} '_N']) = zeros(1, numTop);
end

for i=1:length(testInd)
    curLbl = currentLabel{i};
    fprintf('# Original: %s\n', cell2mat(currentLabel{i}));
    fprintf('# Predicted: X  -  Y:\n');
    if testOnTraining
        % Remove 1st result, which is (normally) the one already given
        fprintf('  %s - %s \n\n',cell2mat(indicesX{i}(3:end)),cell2mat(Results{i}(3:end)))
        predLbls = Results{i}(3:2:3+numTop);
    else
        fprintf('  %s - %s \n\n',cell2mat(indicesX{i}),cell2mat(Results{i}))
        predLbls = Results{i}(1:2:1+numTop);
    end
    curLbl = curLbl{1}(1); 
    for lb = 1:length(predLbls)
        predLbls{lb} = predLbls{lb}(1);
    end
    if strcmp(curLbl, 't'), curLbl = 'o'; end % treat objects and tools as one class
       
    accuracyMRD.(curLbl) = accuracyMRD.(curLbl) + sum(strcmp(curLbl, predLbls).*weights);
    accuracyMRD.([curLbl '_N']) = accuracyMRD.([curLbl '_N']) + 1;
    
    randVec = rand(1,numTop);
    for lb = 1:numTop
        if randVec(lb) <= sum(priorOnClasses(1))
            predRand{lb} = classes{1};
        elseif randVec(lb) <= sum(priorOnClasses(1:2))
            predRand{lb} = classes{2};
        elseif randVec(lb) <= sum(priorOnClasses(1:3))
            predRand{lb} = classes{3};
        end
    end
    
    accuracyRAND.(curLbl) = accuracyRAND.(curLbl) + sum(strcmp(curLbl, predRand)).*weights;
    accuracyRAND.([curLbl '_N']) = accuracyRAND.([curLbl '_N']) + 1;
end

%%
for i=1:length(classes)
    curLbl = classes{i};
    accuracyMRD.(curLbl)  = accuracyMRD.(curLbl)  / accuracyMRD.([curLbl '_N']) ;
    accuracyRAND.(curLbl)  = accuracyRAND.(curLbl) / accuracyRAND.([curLbl '_N']) ;
end
%% %%%

%{ 
% OLD----
if testOnTraining
    % Find p(X_* | Y_*): If Y* is not taken from the tr. data, then this step
    % must be an optimisation step of q(x*). X_* and Y_* here refer to the
    % spaces of the submodel obsMod.
    y_star = model.comp{obsMod}.y(i);
    x_star = model.comp{obsMod}.vardist.means(i,:);
    varx_star = model.comp{obsMod}.vardist.covars(i,:);
    % Find the 5 closest latent points to the x_*, based only on the
    % sharedDimensions. Since we test on a training point, one of the
    % distances is going to be 0, i.e. the test point itself. What is
    % interesting to see, is the 2nd,3rd .... ranked x's, they should be
    % different but showing that the sharedDims are the ones that model the
    % ambiguities, i.e. the different animals perhaps (i.e we would like to
    % see only animals or only objects).
    [ind, distInd] = nn_class(model.X(:,sharedDims), x_star(:,sharedDims), 5, 'euclidean');
    
    ZpredMu = zeros(length(ind), size(model.comp{infMod}.y,2));
    ZpredSigma = zeros(length(ind), size(model.comp{infMod}.y,2));
    
    % Find p(y_*|x_*) for every x_* found from the NN
    for k=1:length(ind)
        [ZpredMu(k,:), ZpredSigma(k,:)] = vargplvmPosteriorMeanVar(model.comp{infMod}, model.X(ind(k),:));
    end
    
    indYnn = [];
    % Do a NN on the DATA space, for every predicted output.
    for j=1:size(ZpredMu,1)
        [indYnn(j), distInd] = nn_class(model.comp{infMod}.y, ZpredMu(j,:),1,'euclidean');
    end
    for j=1:length(indYnn)
        activity{infMod}(indYnn(j))
    end
    'Original:'
    activity{obsMod}(i,:)
else
    testInd = 5:30;
    Yts = YallTest{obsMod}(testInd,:); % Now this is a matrix
    for i=1:size(Yts,1)
        %{
        % initialize the latent points using the nearest neighbour
        % from the training data
        dst = dist2(Yts(i,:), model.comp{obsMod}.y);
        [mind, mini] = min(dst);
        fprintf('# Initialise using label %s\n', cell2mat(activity{obsMod}(mini)));
        %}
        currentLabel = activityTest{obsMod}(testInd(i));
        fprintf('# Original: %s\n', cell2mat(currentLabel));
        % Initialize the latent point using a point of the training data which has the same label
        mini = find(strcmp(activity{obsMod},currentLabel));
        mini = mini(1);
        Init(i,:) = model.vardist.means(mini,:);
        fprintf('# Init. vardist with y with label: %s\n', cell2mat(activity{obsMod}(mini)));
        vardistx = vardistCreate(model.comp{obsMod}.vardist.means(mini,:), model.q, 'gaussian');
        vardistx.covars = model.comp{obsMod}.vardist.covars(mini,:);
        model.comp{obsMod}.vardistx = vardistx;
        display=1;
        iters = 3;
        % Find p(X_* | Y_*) which is approximated by q(X_*)
        [x_star, varx_star, modelUpdated] = vargplvmOptimisePoint(model.comp{obsMod}, vardistx, Yts(i,:), display, iters);%%%
        numberOfNN = 5;
        % Now we optimised a datapoint X_* by taking into account all
        % dimensions. Now, based on the shared dimensions of
        % that, we select the closest (in a NN manner) X from the training
        % data.
        [ind, distInd] = nn_class(model.X(:,sharedDims), x_star(:,sharedDims), 5, 'euclidean');
        
        ZpredMu = zeros(length(ind), size(model.comp{infMod}.y,2));
        ZpredSigma = zeros(length(ind), size(model.comp{infMod}.y,2));
        fprintf('# NN labels for X:\n');
        % Find p(y_*|x_*) for every x_*
        % Careful:  Since x is coming from the training set, it is somewhat
        % natural that the predicted y's will also closely match the ones
        % from the training set. What makes it interesting, is to show that
        % the shared dimensions are the ones controlling the common
        % ambiguities, so, actually it really matters what indices of
        % training X's we find with the NN based on sharedDims.
        for k=1:length(ind)
            [ZpredMu(k,:), ZpredSigma(k,:)] = vargplvmPosteriorMeanVar(model.comp{infMod}, model.X(ind(k),:));
            fprintf('   %s\n',cell2mat(activity{obsMod}(ind(k))));
        end
        
        indYnn = [];
        % Do a NN on the DATA space, for every predicted output.
        for j=1:size(ZpredMu,1)
            [indYnn(j), distInd] = nn_class(model.comp{infMod}.y, ZpredMu(j,:),1,'euclidean');
        end
        
        fprintf('# Predicted:\n');
        Results{i} = [];
        for j=1:length(indYnn)
            Results{i} = [Results{i} ',' activity{infMod}(indYnn(j))];
            fprintf('  %s\n',cell2mat(activity{infMod}(indYnn(j))));
        end
    end
    
    % Now compare activityTest{obsMod}(testInd(i)) with Results{i} for
    % every i
    for i=1:size(Yts,1)
        fprintf('# Original: %s\n', cell2mat(activityTest{obsMod}(testInd(i))));
        fprintf('# Predicted: ');
        fprintf('%s\n\n',cell2mat(Results{i}(2:end)))
    end
    %errsumFull = sum((ZpredMu - Z_test(testInd,:)).^2);
    %errorFull = mean(errsumFull);
end
%}