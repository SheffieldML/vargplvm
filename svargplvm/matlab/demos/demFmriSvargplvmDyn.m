% DEMFMRISVARGPLVMDYN Run the shared variational GPLVM on fmri data by also
% constraining with labels
% 
% COPYRIGHT: Andreas C. Damianou
%
% SHEFFIELDML


% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

% Define constants (in a manner that allows other scripts to parametrize
% this one).
if ~exist('experimentNo') ,  experimentNo = 300 ;      end
if ~exist('itNo')         ,  itNo = [700 1500 500 500];              end     % Default: 2000
if ~exist('indPoints')    ,  indPoints = 80;          end     % Default: 49
if ~exist('latentDims')    ,  latentDim = 12;          end
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
% Options: 'concatenated', 'separately'
if ~exist('initial_X'), initial_X = 'concatenated'; end 
if ~exist('dataType'), dataType = 'fmriWithLbls'; end
if ~exist('enableParallelism'), enableParallelism = 1; end
if ~exist('doPredictions'), doPredictions = 0; end
if ~exist('keepDiary'), keepDiary = false; end
%%
dynamicsConstrainType = {'labels'};
DgtN = 1;


% This script initialises the options structure 'globalOpt'.
svargplvm_init;

globalOpt.dataSetNames = {'fmri400', 'fmri900'};
globalOpt.dataSetName = 'fmri';


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
lblsTrain = zeros(N,numClasses);
% 1-of K encoding: animal is 100, other is 010, rest is 001
lblsTrain(animalIndTr,1) = 1;
lblsTrain(otherIndTr,2) = 1;
lblsTrain(restIndTr,3) = 1;
labelsTrain = transformLabels(lblsTrain);

%% 


options = svargplvmOptions(Ytr, globalOpt, labelsTrain);

if ~isempty(globalOpt.dynamicsConstrainType)
    for i=1:numberOfDatasets
        % Set up dynamcis (i.e. back-constraints) model
        optionsDyn{i}.type = 'vargpTime';
        optionsDyn{i}.inverseWidth=30;
        %   optionsDyn.vardistCovars = vardistCovarsMult;
        optionsDyn{i}.initX = globalOpt.initX;
        optionsDyn{i}.constrainType = globalOpt.dynamicsConstrainType;
        
        if exist('timeStampsTraining')
            optionsDyn{i}.t = timeStampsTraining;
        end
        if exist('labelsTrain') && ~isempty(labelsTrain)
            optionsDyn{i}.labels = labelsTrain;
        end
    end
else
    optionsDyn= [];
end

%%
model = svargplvmModelCreate(Ytr, globalOpt, options, optionsDyn);
model.globalOpt = globalOpt;
model.options = options;
       
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
    try
         pool_open = matlabpool('size')>0;
    catch e
        pool_open = 0;
    end
    if ~pool_open
        warning('enableParallelism active but matlabpool is closed!')
    end
    
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
%%
model = svargplvmOptimiseModel(model);

svargplvmShowScales(model)







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

Nstar = 100;
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
