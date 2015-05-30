% DEMYALESVARGPLVM4 Run the Shared Var. GP-LVM on a subset of the Yale
% faces.
% DESC Run the Shared Var. GP-LVM on a subset of the Yale faces. The code
% for creating this subset out of raw images exists in comments. Unlike
% demYaleSvargplvm1, this demo is not a wrapper, it can be used as a
% standalone demo.
%
% VARGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);


% Create the dataset out of the images.
baseDir=[localDatasetsDirectorySmall 'sharedVargplvm' filesep 'artificialLight'];
selDirs = {'y1','y2','y3','y4','y5','y6'};

for d=1:length(selDirs)
    dirFrom=[baseDir filesep selDirs{d}];
    a=dir(dirFrom);
    counter = 0;
    for i=1:length(a)
        if length(a(i).name)>4 & strcmp(a(i).name(end-2:end),'pgm')
            im = imread([dirFrom filesep a(i).name]);
            %imagesc(im), colormap('gray'); title(a(i).name), pause
            counter = counter+1;
            Yall{d}(counter,:)=im(:)';
        end
    end
    Yall{d} = double(Yall{d});
end
height = size(im,1);
width = size(im,2);
for d=1:length(Yall)
    Yall{d} = util_artificialShadow(reshape(Yall{d},height,width),6,30,[],3000);
end


%{
for i=1:size(Yall{1},1)
    for d=1:length(Yall)
        subplot(2,3,d)
        imagesc(reshape(Yall{d}(i,:), height, width)); colormap('gray');
    end
    if i==1
        pause
    else
        pause(0.1)
    end
end
%}

Y{1} = Yall{1};
Y{2} = Yall{2};
Y{3} = Yall{3};
Y{4} = [];
Y{5} = [];
Y{6} = [];
identities{1}=[ones(size(Y{1},1),1); 2*ones(size(Y{2},1),1); 3*ones(size(Y{3},1),1)];
idTmp{4}=[];idTmp{5}=[];idTmp{6}=[];
for i=1:size(Yall{4},1)
    perm = randperm(3)+3;%perm=[1 2 3]+3;
    Y{4} = [Y{4}; Yall{perm(1)}(i,:)]; idTmp{4}=[idTmp{4}; perm(1)];
    Y{5} = [Y{5}; Yall{perm(2)}(i,:)]; idTmp{5}=[idTmp{5}; perm(2)];
    Y{6} = [Y{6}; Yall{perm(3)}(i,:)]; idTmp{6}=[idTmp{6}; perm(3)];
end
clear 'Yall'
Yall{1} = [Y{1}; Y{2}; Y{3}];
Yall{2} = [Y{4}; Y{5}; Y{6}];
identities{2} = [idTmp{4};idTmp{5};idTmp{6}];
clear 'Y' 'idTmp' 'd'

numberOfDatasets = length(Yall);

%{
for i=1:size(Yall{1},1)
    for d=1:length(Yall)
        subplot(1,2,d)
        imagesc(reshape(Yall{d}(i,:), height, width)); colormap('gray'); title(num2str(identities{d}(i)))
    end
    if i==1
        pause
    else
        pause%(0.3)
    end
end
%}
%%

if ~exist('trainModel','var'), trainModel = true; end
if ~exist('itNo')         ,  itNo = [500 1500 1500];              end     % Default: 2000
if ~exist('indPoints')    ,  indPoints = 45;          end     % Default: 49
if ~exist('initVardistIters'), initVardistIters = 800;      end
if ~exist('mappingKern')   ,  mappingKern = 'rbfardjit'; end

% Set to empty value {} to work with toy data
if ~exist('dataSetNames')    ,    dataSetNames = {};    end
if ~exist('invWidthMult'),       invWidthMult = 5;                     end
if ~exist('dataType'), dataType = 'artificialLight'; end
if ~exist('latentDimPerModel'), latentDimPerModel = 7; end
if ~exist('experimentNo'), experimentNo = 1; end
if ~exist('doPredictions'), doPredictions = false; end
% If this is true, then the model is in "D > N" mode.
if ~exist('DgtN'), DgtN = true; end
% Create initial X by doing e.g. ppca in the concatenated model.m's or by
% doing ppca in the model.m's separately and concatenate afterwards?
if ~exist('initial_X'), initial_X = 'separately'; end % Other options: 'concatenated'
% Which indices to use for training, rest for test
if ~exist('indTr'), indTr = -1; end

enableParallelism = 0;

if exist('pyramid','var') % extract pyramid representation of the images
    if pyramid
        for e=1:size(Y,2)
            Y{e} = im2pyramid(Y{e}, lbls(1), lbls(2), 4);
        end
    end
end





%-- Load datasets
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
if trainModel
    
    %-- Options for the models
    for i=1:numberOfDatasets
        % Set up models
        options{i} = vargplvmOptions('dtcvar');
        options{i}.kern = mappingKern; %{'rbfard2', 'bias', 'white'};
        %indPoints = 80; %%%%%
        options{i}.numActive = indPoints;
        options{i}.optimiser = 'scg2';
        if ~DgtN
            options{i}.enableDgtN = false;
        end
        % !!!!! Be careful to use the same type of scaling and bias for all
        % models!!!
        
        % scale = std(Ytr);
        % scale(find(scale==0)) = 1;
        %options.scaleVal = mean(std(Ytr));
        options{i}.scaleVal = sqrt(var(Ytr{i}(:)));
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
        
        %    mAll = [mAll m{i}]; % Concatenation (doesn't work if different sizes)
    end
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
        
        
        
        model{i}.beta=1/(0.01*var(m{i}(:)));
        prunedModelInit{i} = vargplvmPruneModel(model{i});
        %disp(model{i}.vardist.covars)
    end
    
    
    
    %modelInit = model;%%%TEMP
    
    %--  Unify models into a structure
    svargplvm_init
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
    
    %%
    
    
    display = 1;
    %%%% Optimisation
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
else
    capName = dataType;
    capName(1) = upper(capName(1));
    modelType = 'svargplvm';
    modelType(1) = upper(modelType(1));
    fileToSave = ['dem' capName modelType num2str(experimentNo) '.mat'];
    load(fileToSave);
    model = svargplvmRestorePrunedModel(prunedModel, Ytr);clear('prunedModel');clear('prunedModelInit');
end

%%

if ~exist('resultsDynamic')  resultsDynamic = 0; end


for i=1:length(model.comp)
    if model.comp{i}.DgtN
        model.comp{i}.m = model.comp{i}.mOrig;
    end
end
for i=1:length(model.comp)
    model.comp{i}.vis.startDim = {1,2};
    model.comp{i}.vis.startPos = model.vardist.means(1,:);
end

%%
v = 1;
modelVis = model.comp{v};
if resultsDynamic
    %bar(model.comp{v}.kern.comp{1}.inputScales);
    figure
    % The following causes OUTOFMEMORY exception except from when we prune the
    % video dimensions: (Note: also, lvmVisualise was changed a bit so that no
    % dynamic slides is presented, because otherwise a strange error occurs).
    modelVis.y = Ytr{v};
    [modelP, newHeight, newWidth] = vargplvmReduceVidModel(modelVis, height, width, 2,2);
    lvmVisualiseGeneral(modelP, [], 'imageMRDVisualise', 'imageMRDModify', false, [newHeight newWidth],0,0,1);
    clear modelP
end


%%

%---------------------------- PREDICTIONS ---------------
if ~doPredictions
    return
end


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

if ~exist('testOnTraining')
    testOnTraining=1;
end

numberTestPoints = 10;
if testOnTraining
    perm = randperm(model.N);
    testInd = perm(1:numberTestPoints);
else
    perm = randperm(size(Yts{obsMod},1));
    testInd = perm(1:numberTestPoints);
end

scrsz = get(0,'ScreenSize');

for i=1:length(testInd)
    curInd = testInd(i);
    fprintf('# Testing indice number %d ', curInd);
    if testOnTraining
        fprintf('taken from the training set\n');
        y_star = model.comp{obsMod}.y(curInd,:);
        x_star = model.comp{obsMod}.vardist.means(curInd,:);
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
        [x_star, varx_star, modelUpdated] = vargplvmOptimisePoint(model.comp{obsMod}, vardistx, y_star, display, iters);%%%
    end
    numberOfNN = 9;
    % Now we selected a datapoint X_* by taking into account only the
    % private dimensions for Y. Now, based on the shared dimensions of
    % that, we select the closest (in a NN manner) X from the training data.
    fprintf('# Finding the %d NN of X_* with the training X based only on the shared dims.\n', numberOfNN);
    [ind, distInd] = nn_class(model.X(:,sharedDims), x_star(:,sharedDims), numberOfNN, 'euclidean');
    
    ZpredMu = zeros(length(ind), size(model.comp{infMod}.y,2));
    ZpredSigma = zeros(length(ind), size(model.comp{infMod}.y,2));
    
    
    % Find p(y_*|x_*) for every x_* found from the NN
    fprintf('# Predicting images from the NN of X_* ');
    for k=1:numberOfNN
        fprintf('.');
        x_cur = model.X(ind(k),:);
        %x_cur(sharedDims) = x_star(sharedDims); %%% OPTIONAL!!!
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
    
    
    %{
    indYnn = [];
    % Do a NN on the DATA space, for every predicted output.
    for j=1:size(ZpredMu,1)
        [indYnn(j), distInd] = nn_class(model.comp{infMod}.y, ZpredMu(j,:),1,'euclidean');
    end
    for j=1:length(indYnn)
        figure, imagesc(reshape(model.comp{infMod}.y(indYnn(j),:),height, width)), title([num2str(j)]), colormap('gray')
    end
    %}
end

if ~testOnTraining
    errsumFull = sum((ZpredMu - Yts).^2);
    errorFull = mean(errsumFull);
end


%{
% OLD CODE!!
if testOnTraining
  %  i = size(model.comp{obsMod}.y,1); % last observation
   i=23;
   %i=120;
    % Find p(X_* | Y_*): If Y* is not taken from the tr. data, then this step
    % must be an optimisation step of q(x*). X_* and Y_* here refer to the
    % spaces of the submodel obsMod.
    y_star = model.comp{obsMod}.y(i);
    x_star = model.comp{obsMod}.vardist.means(i,:);
    varx_star = model.comp{obsMod}.vardist.covars(i,:);
    % Find the 10 closest latent points to the x_*, based only on the
    % sharedDimensions (since we test on a training point, one of the
    % distances is going to be 0, i.e. the test point itself).
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
        figure, imagesc(reshape(model.comp{infMod}.y(indYnn(j),:),height, width)), title([num2str(j)]), colormap('gray')
    end
    figure, imagesc(reshape(model.comp{obsMod}.y(i,:),height,width)), title('Original'), colormap('gray')
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
        mu = vargplvmPosteriorMeanVar(model.comp{infMod}, x);
        
        ZpredK{i} = mu;
        
        ZpredMu(i,:) = mu(1);
        
        %ZpredSigma(i,:) = sigma;
    end
    
    errsumFull = sum((ZpredMu - Z_test(testInd,:)).^2);
    errorFull = mean(errsumFull);
end

%}