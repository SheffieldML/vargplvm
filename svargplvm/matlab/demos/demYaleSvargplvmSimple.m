% DEMYALESVARGPLVMSIMPLE Run the Shared Var. GP-LVM on a subset of the Yale
% faces.
% DESC Run the Shared Var. GP-LVM on a subset of the Yale faces. This demo
% is similar to demYaleVargplvm2.m but with fewer options.

% COPYRIGHT: Andreas C. Damianou, Carl Henrik Ek,  2011
% SEEALSO : demYaleSvargplvm2

% SHEFFIELDML


% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

%------------ Constants -------------
if ~exist('itNo')         ,  itNo = [500 1500 1500 1500];     end     % Default: 2000
if ~exist('indPoints')    ,  indPoints = 120;          end     % Default: 49
if ~exist('initVardistIters'), initVardistIters = 180;      end
if ~exist('mappingKern')   ,  mappingKern = 'rbfardjit'; end
if ~exist('latentDimPerModel'), latentDimPerModel = 10; end
if ~exist('experimentNo'), experimentNo = 404; end
% Set to true to explore the latent space after optimising.
if ~exist('resultsDynamic')  resultsDynamic = 0; end
% Set to true to retrain the model, otherwise a already trained model is
% loaded.
reTrainModel = false;
% Indices (datapoints) to be used for training (the rest form a test set).
indTr = [1:3:39 40:192];
dataType = 'Yale6Sets';
dataSetNames = 'YaleSubset6_1';
%-------------------------------------

[Y,lbls]=svargplvmLoadData(dataSetNames);

%------- Shuffle data ---
N1 = size(Y{1},1);
Yall{2} = [Y{4};Y{5};Y{6}];
identities{2}=[ones(N1,1) 2*ones(N1,1) 3*ones(N1,1)];

numSubsets = 3;
Yall{1} = zeros(numSubsets*N1, size(Y{1},2));
for i=1:N1
    perm = randperm(numSubsets);
    counter = 0;
    for j=perm
        Yall{1}(i+counter*N1,:) = Y{j}(i,:);
        identities{1}(i+counter*N1) = j;
        counter = counter+1;
    end
end
%------------------------

clear Y;


numberOfDatasets = length(Yall);
height = lbls(1); width = lbls(2);


%{
%-- Visualise the dataset
for i=1:size(Yall{1},1)
    for d=1:numberOfDatasets
        subplot(1,numberOfDatasets,d)
        imagesc(reshape(Yall{d}(i,:),height, width)), title(num2str(identities{d}(i))),colormap('gray');
    end
    pause
end
%}




%-- Reform datasets and split into training and test set
for i=1:numberOfDatasets
    Y = Yall{i};
    dims{i} = size(Y,2);
    N{i} = size(Y,1);
    if indTr == -1
        indTr = 1:N{i};
    end
    indTs = setdiff(1:size(Y,1), indTr);
    Ytr{i} = Y(indTr,:);
    Yts{i} = Y(indTs,:);
    t{i} = linspace(0, 2*pi, size(Y, 1)+1)'; t{i} = t{i}(1:end-1, 1);
    timeStampsTraining{i} = t{i}(indTr,1); %timeStampsTest = t(indTs,1);
    d{i} = size(Ytr{i}, 2);
end

for i=2:numberOfDatasets
    if N{i} ~= N{i-1}
        error('The number of observations in each dataset must be the same!');
    end
end

%%
if ~reTrainModel
    load demYale6SetsSvargplvm4
    model = svargplvmRestorePrunedModel(prunedModel, Ytr);
else
    %-- Options for the models
    for i=1:numberOfDatasets
        % Set up models
        options{i} = vargplvmOptions('dtcvar');
        options{i}.kern = mappingKern; %{'rbfard2', 'bias', 'white'};
        options{i}.numActive = indPoints;
        options{i}.optimiser = 'scg2';
        options{i}.scaleVal = sqrt(var(Ytr{i}(:)));
        options{i}.enableDgtN = false;

    end
    
    %-------------- INIT LATENT SPACE ---%
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
    
    fprintf('# Initialising X by performing ppca in each observed (scaled) dataset separately and then concatenating...\n');
    X_init{1} = ppcaEmbed(m{1},latentDimPerModel);
    X_init{2} = ppcaEmbed(m{2},latentDimPerModel);
    X_init = [X_init{1} X_init{2}];
    
    %-----------------
    
    latentDim = size(X_init,2);
    
    % Free up some memory
    clear('Y')
    
    
    
    %-- Create the sub-models.
    for i=1:numberOfDatasets
        fprintf(1,'# Creating the model...\n');
        options{i}.initX = X_init;
        model{i} = vargplvmCreate(latentDim, d{i}, Ytr{i}, options{i});
        model{i}.X = X_init;
        model{i} = vargplvmParamInit(model{i}, model{i}.m, model{i}.X);
        model{i}.X = X_init;
        
        inpScales = 5./(((max(model{i}.X)-min(model{i}.X))).^2); % Default 5
        model{i}.kern.comp{1}.inputScales = inpScales;
        
        if strcmp(model{i}.kern.type, 'rbfardjit')
            model{i}.kern.inputScales = model{i}.kern.comp{1}.inputScales;
        end
        params = vargplvmExtractParam(model{i});
        model{i} = vargplvmExpandParam(model{i}, params);
        model{i}.vardist.covars = 0.5*ones(size(model{i}.vardist.covars)) + 0.001*randn(size(model{i}.vardist.covars));
        
        
        
        model{i}.beta=1/(0.01*var(model{i}.m(:)));
        prunedModelInit{i} = vargplvmPruneModel(model{i});
        %disp(model{i}.vardist.covars)
    end
    
    
    
    %--  Unify models into a structure
    model = svargplvmModelCreate(model);
    model.dataSetNames = dataSetNames;
    model.experimentNo = experimentNo;
    model.dataType = dataType;
    %%---
    
    capName = dataType;
    capName(1) = upper(capName(1));
    modelType = model.type;
    modelType(1) = upper(modelType(1));
    fileToSave = ['dem' capName modelType num2str(experimentNo) '.mat'];
    %%---
    
    
    % Force kernel computations
    params = svargplvmExtractParam(model);
    model = svargplvmExpandParam(model, params);
    
    %%
    
    %-------------  Optimisation ---------------
    display = 1;
    % do not learn beta and sigma_f for few iterations for intitialization
    if initVardistIters ~=0
        model.initVardist = 1; model.learnSigmaf = 0;
        model = svargplvmPropagateField(model,'initVardist', model.initVardist);
        model = svargplvmPropagateField(model,'learnSigmaf', model.learnSigmaf);
        fprintf(1,'# Intitiliazing the variational distribution...\n');
        model = svargplvmOptimise(model, display, initVardistIters); % Default: 20
        %fprintf(1,'1/b = %.4d\n',1/model.beta);
        
        modelInitVardist = model;
        model.initVardistIters=initVardistIters;
    end
    
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
end
%%
%------------- Explore latent space

v = 1; % This can be either 1 or 2
modelVis = model.comp{v};
if resultsDynamic
    %bar(model.comp{v}.kern.comp{1}.inputScales);
    figure
    % The following causes OUTOFMEMORY exception except from when we prune the
    % video dimensions: (Note: also, lvmVisualise was changed a bit so that no
    % dynamic slides is presented, because otherwise a strange error occurs).
    modelVis.y = Ytr{v};
    sc = 2; % Set to 4 and try again, if you run out of memory
    [modelP, newHeight, newWidth] = vargplvmReduceVidModel(modelVis, height, width, sc,sc);
    lvmVisualiseGeneral(modelP, [], 'imageMRDVisualise', 'imageMRDModify', false, [newHeight newWidth],0,0,1);
    clear modelP
end

%%






%%
doPredictions = 1;
%---------------------------- PREDICTIONS ---------------
if ~doPredictions
    return
end

% Set to 1 to test on the training data itself, set to 0 to use the test
% dataset.
if ~exist('testOnTraining')
    testOnTraining=1;
end


obsMod = 1; % one of the involved sub-models (possible values: 1 or 2).
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


% Number of test points to use
numberTestPoints = 15;
if testOnTraining
    perm = randperm(model.N);
    testInd = perm(1:numberTestPoints);
else
%    perm = randperm(size(Yts{obsMod},1));
%    testInd = perm(1:numberTestPoints);
    testInd = 1:size(Yts{1},1);
end

scrsz = get(0,'ScreenSize');

x_star = zeros(length(testInd), size(model.X,2));
for i=1:length(testInd)
    curInd = testInd(i);
    fprintf('# Testing indice number %d ', curInd);
    if testOnTraining
        fprintf('taken from the training set\n');
        y_star = model.comp{obsMod}.y(curInd,:);
        x_star(i,:) = model.comp{obsMod}.vardist.means(curInd,:);
        varx_star = model.comp{obsMod}.vardist.covars(curInd,:);
    else
        fprintf('taken from the test set\n');
        y_star = Yts{obsMod}(curInd,:);
        z_star = Yts{infMod}(curInd,:);
        dst = dist2(y_star, model.comp{obsMod}.y);
        [mind, mini] = min(dst);
        
        Init(i,:) = model.vardist.means(mini,:);
        vardistx = vardistCreate(model.comp{obsMod}.vardist.means(mini,:), model.q, 'gaussian');
        vardistx.covars = model.comp{obsMod}.vardist.covars(mini,:);
        model.comp{obsMod}.vardistx = vardistx;
        display=1;
        iters = 250;
        % Find p(X_* | Y_*) which is approximated by q(X_*)
        [x_star(i,:), varx_star, modelUpdated] = vargplvmOptimisePoint(model.comp{obsMod}, vardistx, y_star, display, iters);%%%
    end
    numberOfNN = 9;
    % Now we selected a datapoint X_* by taking into account only the
    % private dimensions for Y. Now, based on the shared dimensions of
    % that, we select the closest (in a NN manner) X from the training data.
    fprintf('# Finding the %d NN of X_* with the training X based only on the shared dims.\n', numberOfNN);
    
    % Actually this is a weighted NN.
    w = s1(sharedDims)+s2(sharedDims);
    [ind, distInd] = nn_class2(model.X(:,sharedDims), x_star(i,sharedDims), numberOfNN, 'weighted', w);
    %[ind, distInd] = nn_class(model.X(:,sharedDims), x_star(i,sharedDims), numberOfNN, 'euclidean');
    
    ZpredMu = zeros(length(ind), size(model.comp{infMod}.y,2));
    ZpredSigma = zeros(length(ind), size(model.comp{infMod}.y,2));
    
    
    % Find p(y_*|x_*) for every x_* found from the NN
    fprintf('# Predicting images from the NN of X_* ');
    for k=1:numberOfNN
        fprintf('.');
        x_cur = model.X(ind(k),:);
        % Make the shared dimensions of the current NN the same as the
        % ones of the latent point for y_star
       %  x_cur(sharedDims) = x_star(i,sharedDims); %%% OPTIONAL!!!
        
        % Make the shared dimensions for the current NN the same as the 
        % ones of the closest NN to y_star but from the Z dataset
       % x_cur(sharedDims) = model.X(ind(1),sharedDims); %%%%% OPTIONAL #2 !!!!
        
        %[ZpredMu(k,:), ZpredSigma(k,:)] = vargplvmPosteriorMeanVar(model.comp{infMod}, model.X(ind(k),:));
        ZpredMu(k,:) = vargplvmPosteriorMeanVar(model.comp{infMod}, x_cur);
    end
    fprintf('\n\n');
    
    
    %-- Plots

    % Open a big figure (first 2 args control the position, last 2 control
    % the size)
    figure('Position',[scrsz(3)/100.86 scrsz(4)/6.666 scrsz(3)/1.0457 scrsz(4)/1.0682],...
        'Name',['Fig: ' num2str(i) ' (Exp: ' num2str(experimentNo) ')'],'NumberTitle','off')
    numRows = 2;
    
    if testOnTraining
        numCols = ceil((numberOfNN+1)/numRows);
        plotCounter = 1;
    else
        % For the real test image!
        numCols = ceil((numberOfNN+2)/numRows);
        plotCounter = 2;
    end
    subplot(numRows, numCols, 1)
    imagesc(reshape(y_star,height,width)), title(['Original y (image #' num2str(curInd) ')']), colormap('gray')
    
    if ~testOnTraining
        subplot(numRows, numCols, 2)
        imagesc(reshape(z_star,height,width)), title(['Corresponding z (image #' num2str(curInd) ')']), colormap('gray')
    end
    
    for k=1:numberOfNN
        subplot(numRows, numCols, k+plotCounter)
        imagesc(reshape(ZpredMu(k,:), height, width)), title(['NN #' num2str(k)]), colormap('gray');
    end
end


