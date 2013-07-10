% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

% Define constants (in a manner that allows other scripts to parametrize
% this one).
if ~exist('experimentNo') ,  experimentNo = 404;      end
if ~exist('itNo')         ,  itNo = [600 600 600];   end     % Default: 2000
if ~exist('indPoints')    ,  indPoints = 64;          end     % Default: 49
if ~exist('latentDim')    ,  latentDim = 7;          end
% Set to 1 to use dynamics or to 0 to use the standard var-GPLVM
if ~exist('dynUsed')      ,  dynUsed = 0;             end
if ~exist('initVardistIters'), initVardistIters = 180;      end     % DEFAULT: 23
if ~exist('dynamicKern')   ,     dynamicKern = {'rbf', 'white', 'bias'}; end
if ~exist('mappingKern')   ,     mappingKern = 'rbfardjit'; end %{'rbfard2', 'bias', 'white'}; end
if ~exist('vardistCovarsMult'),  vardistCovarsMult=1.3;                  end
if ~exist('invWidthMultDyn'),    invWidthMultDyn = 100;                     end
if ~exist('invWidthMult'),       invWidthMult = 5;                     end
if ~exist('initX'),     initX ='ppca';   end
if ~exist('DgtN'), DgtN = 1; end

svargplvmModel = 'demYale6SetsSvargplvm25';
if ~exist('svargplvmYtr')
    %svargplvmYtr = load(''); % ASSUME it's already loaded.
    load svargplvmYtr
end

% Create the dataset out of the images.
baseDir=[localDatasetsDirectoryLarge 'CroppedYale' filesep 'CroppedYale'];
if ~exist('selDir')
    selDir = '03';  %{'03', '08'}; % 03,08
end

dirFrom=[baseDir filesep 'yaleB' selDir];
a=dir(dirFrom);
counter = 0;
for i=1:length(a)
    if length(a(i).name)>4 & strcmp(a(i).name(end-2:end),'pgm') ...
            & ~strcmp(a(i).name(end-10:end-4),'Ambient')
        im = imread([dirFrom filesep a(i).name]);
        %imagesc(im), colormap('gray'); title(a(i).name), pause
        counter = counter+1;
        Y(counter,:)=im(:)';
    end
end
Y = double(Y);
height = size(im,1);
width = size(im,2);
clear('im','a');

%-- Split into training and test set
if ~exist('indTs')
   % indTs = [17,22,27,39,45,49]; %[5,11,14];
   indTs = [];
end
indTr = setdiff( 1:size(Y,1),indTs);
Ytr = Y(indTr,:);
Yts = Y(indTs, :);
indPoints = length(indTr);
%---
d = size(Ytr,2);

dataSetName = ['YaleFace' selDir];
%{
    [Y, lbls] = vargplvmLoadData(dataSetName);
    height = lbls(1); width = lbls(2);

    %%
    for i=1:size(Y,1)
        imagesc(reshape(Y(i,:),height, width)), title(num2str(i)),colormap('gray');
        pause
    end
%}

%%

%----
% Load svargplvm model and find the common latent space.
load(svargplvmModel);
modelS = prunedModel;
if ~exist('sharedDims')
    s1 = modelS.comp{1}.kern.comp{1}.inputScales;
    s2 = modelS.comp{2}.kern.comp{1}.inputScales;
    % Normalise values between 0 and 1
    s1 = s1 / max(s1);
    s2 = s2 / max(s2);  
    %  thresh = max(model.comp{obsMod}.kern.comp{1}.inputScales) * 0.001;
    thresh = 0.005;
    retainedScales{1} = find(s1 > thresh);
    %thresh = max(model.comp{infMod}.kern.comp{1}.inputScales) * 0.001;
    retainedScales{2} = find(s2  > thresh);
    sharedDims = intersect(retainedScales{1}, retainedScales{2});
end
sharedDims = [1 2 3];%%%%%%%%
svargplvmX = modelS.X;
svargplvmVardist = modelS.vardist;
clear('modelS','prunedModel','prunedModelInit')
%----

% Set up model
options = vargplvmOptions('dtcvar');
options.kern = mappingKern; %{'rbfard2', 'bias', 'white'};
options.numActive = indPoints;
options.optimiser = 'scg2'; % Set 'scg' if 'scg2' is missing
if ~exist('DgtN') || ~DgtN
    options.enableDgtN = false;
end
options.scaleVal = sqrt(var(Ytr(:)));


% Compute m, the normalised version of Ytr (to be used for
% initialisation of X)
bias = mean(Ytr);
scale = ones(1, d);
if(isfield(options,'scale2var1'))
    if(options.scale2var1)
        scale = std(Ytr);
        scale(find(scale==0)) = 1;
        if(isfield(options, 'scaleVal'))
            warning('Both scale2var1 and scaleVal set for GP');
        end
    end
end
if(isfield(options, 'scaleVal'))
    scale = repmat(options.scaleVal, 1, d);
end
% Remove bias and apply scale.
m = Ytr;
for j = 1:d
    m(:, j) = m(:, j) - bias(j);
    if scale(j)
        m(:, j) = m(:, j)/scale(j);
    end
end

% The length(sharedDims) most important dimensions of the inital X for this
% dataset, are going to be replaced by the shared dimensions of the
% svargplvm model's X.
X_init = ppcaEmbed(m,latentDim);


%--- Init. based on the svargplvm mode.
%inpScales = invWidthMult./(((max(X_init)-min(X_init))).^2);
%[void, indScales] = sort(inpScales,'descend');
%tiedX = indScales(1:length(sharedDims));

%tiedX = [1 2 3];%%%%%%%%5

%X_init(:,tiedX) = svargplvmX(1:size(Ytr,1),sharedDims);


options.initX = X_init;

% demo using the variational inference method for the gplvm model
fprintf(1,'# Creating the model...\n');

model = vargplvmCreate(latentDim, d, Ytr, options);
% Temporary: in this demo there should always exist the mOrig field
if ~isfield(model, 'mOrig')
    model.mOrig = model.m;
end
model.X = X_init; %%%%%%%
model = vargplvmParamInit(model, model.mOrig, model.X);
model.X = X_init; %%%%%%%

%-- Make sure that the scales for the tied dimensions are at least as large
% as the largest scale.
model.kern.comp{1}.inputScales = invWidthMult./(((max(model.X)-min(model.X))).^2); % Default 5
%offs = max(model.kern.comp{1}.inputScales) - min(model.kern.comp{1}.inputScales(tiedX)); %%%
%model.kern.comp{1}.inputScales(tiedX) = model.kern.comp{1}.inputScales(tiedX)+offs; %%%%%%%%%
% For the rbfardjit kernel
model.kern.inputScales = model.kern.comp{1}.inputScales;
%--

params = vargplvmExtractParam(model);
model = vargplvmExpandParam(model, params);
%model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));


%--- Init. based on the svargplvm model.
%model.vardist.means(:,tiedX) = svargplvmVardist.means(1:size(Ytr,1), sharedDims);
%model.vardist.covars(:,tiedX) = svargplvmVardist.covars(1:size(Ytr,1), sharedDims);
%model.tiedX = tiedX; %%%%%%%%%


model.beta=1/(0.01*var(model.mOrig(:)));
modelInit = model;

capName = dataSetName;
capName(1) = upper(capName(1));
modelType = model.type;
modelType(1) = upper(modelType(1));
fileToSave = ['dem' capName modelType num2str(experimentNo) '.mat'];
%% ---

display = 1;
%     %%%% Optimisation
%     % do not learn beta for few iterations for intitialization
if initVardistIters ~=0
    model.initVardist = 1;
    model.learnBeta = 0; model.learnSigmaf = 0; % This should be merged with the initVardist field
    fprintf(1,'# Intitiliazing the model (fixed beta)...\n');
    model = vargplvmOptimise(model, display, initVardistIters); % Default: 20
    fprintf(1,'1/b = %.4d\n',1/model.beta);
    model.learnSigmaf = 1; model.learnBeta =1; model.initVardist = 0;
end
model.iters = 0;
prunedModelInit = vargplvmPruneModel(modelInit);
clear modelInit

% Optimise the model.
for i=1:length(itNo)
    iters = itNo(i); % default: 2000
    fprintf(1,'\n# Optimising the model for %d iterations (session %d)...\n',iters,i);
    model = vargplvmOptimise(model, display, iters);
    model.iters = model.iters + iters;
    fprintf(1,'1/b = %.4d\n',1/model.beta);
    modelTr = model;
    fprintf(1,'# 1/b=%.4f\n var(m)=%.4f\n',1/model.beta, var(model.mOrig(:)));
    % Save model
    fprintf(1,'# Saving %s\n',fileToSave);
    prunedModel = vargplvmPruneModel(model);
    save(fileToSave, 'prunedModel', 'prunedModelInit');
end
% Just for compatibility.
if strcmp(model.type,'rbfardjit')
    model.kern.comp{1}.inputScales = model.kern.inputScales;
end

save(fileToSave, 'prunedModel', 'prunedModelInit');

model.m = model.mOrig;
%{
sc = 2; % Set to 4 and try again, if you run out of memory
[modelP, newHeight, newWidth] = vargplvmReduceVidModel(model, height, width, sc,sc);
lvmVisualise(modelP, [], 'imageVisualise', 'imageModify', [newHeight newWidth],0,0,1);
clear modelP
figure,bar(model.kern.comp{1}.inputScales)
%}


%%%%%%%%%TODO ( I can leave 4-5 frames out, and find x* based on y* and
%%%%%%%%%then find the similar pics).


%%% In theory, the sharedDims of the oneFace model might be different. In
%%% this case, I 'll take sharedDimsOneFace = the
%%% length(sharedDimsSvargplvm) largest dims of the vargplvmOne face. Or,
%%% we could do a NN similarity test of all the Xone, Xsvargplvm, to see
%%% the length(sharedDimsSvargplvm) most similar.


doPredictions = 1;
%---------------------------- PREDICTIONS ---------------
if ~doPredictions
    return
end


%%
sharedModel = load(svargplvmModel);
sharedModel = sharedModel.prunedModel;
sharedModel = svargplvmRestorePrunedModel(sharedModel, svargplvmYtr);
sharedModel.comp{1}.m = sharedModel.comp{1}.mOrig;
sharedModel.comp{2}.m = sharedModel.comp{2}.mOrig;

% Set to 1 to test on the training data itself, set to 0 to use the test
% dataset.
if ~exist('testOnTraining')
    testOnTraining=1;
end


% Find the dimensions that are shared for obsMod and infMod
if ~exist('retainedDims')
    snew = model.kern.comp{1}.inputScales;
    % Normalise values between 0 and 1
    snew = snew / max(snew);    
    thresh = 0.08; 
    retainedDims = find(snew > thresh);
end
retainedDims = [2 1 3];%%%%%%%%%%
%retainedDims = [1 2]; sharedDims = [1 2]; %%%%%%%%%%%%%%%%%%%%%%%%%

% Number of test points to use
numberTestPoints = 15;
if testOnTraining
    perm = randperm(model.N);
    testInd = perm(1:numberTestPoints);
else
%    perm = randperm(size(Yts{edobsMod},1));
%    testInd = perm(1:numberTestPoints);
    testInd = 1:size(Yts,1);
end

scrsz = get(0,'ScreenSize');

x_star = zeros(length(testInd), size(model.X,2));
for i=1:length(testInd)
    curInd = testInd(i);
   % fprintf('# Testing indice number %d ', curInd);
    if testOnTraining
    %    fprintf('taken from the training set\n');
        y_star = model.y(curInd,:);
        x_star(i,:) = model.vardist.means(curInd,:);
        varx_star = model.vardist.covars(curInd,:);
    else
     %   fprintf('taken from the test set\n');
        y_star = Yts(curInd,:);
        dst = dist2(y_star, model.y);
        [mind, mini] = min(dst);
        
        Init(i,:) = model.vardist.means(mini,:);
        vardistx = vardistCreate(model.vardist.means(mini,:), model.q, 'gaussian');
        vardistx.covars = model.vardist.covars(mini,:);
        model.vardistx = vardistx;
        display=1;
        iters = 250;
        % Find p(X_* | Y_*) which is approximated by q(X_*)
        [x_star(i,:), varx_star, modelUpdated] = vargplvmOptimisePoint(model, vardistx, y_star, display, iters);%%%
    end
    numberOfNN = 3;
    % Now we selected a datapoint X_* by taking into account only the
    % private dimensions for Y. Now, based on the shared dimensions of
    % that, we select the closest (in a NN manner) X from the training data.
   % fprintf('# Finding the %d NN of X_* with the training X based only on the shared dims.\n', numberOfNN);
    w=s1(sharedDims)+s2(sharedDims);
    [ind, distInd] = nn_class2(sharedModel.X(:,sharedDims), x_star(i,retainedDims), numberOfNN, 'weighted',w);
  % [ind, distInd] = nn_class(sharedModel.X(:,sharedDims), x_star(i,retainedDims), numberOfNN, 'euclidean');
    
    YpredMu = zeros(length(ind), size(sharedModel.comp{1}.y,2));
    ZpredMu = zeros(length(ind), size(sharedModel.comp{2}.y,2));
    
    % Find p(y_*|x_*) for every x_* found from the NN
    fprintf(['# Predicting images from the NN of X_* of ' num2str(curInd) ': ']);
    for k=1:numberOfNN
        curNN = ind(k);
        curNN = curNN - size(Ytr,1)*(ceil(curNN/size(Ytr,1)) -1);
        curNNall(k) = curNN;
        
        fprintf([num2str(curNN) ' ']);
        x_cur = sharedModel.X(ind(k),:);
        % Make the shared dimensions of the current NN the same as the
        % ones of the latent point for y_star
       %  x_cur(sharedDims) = x_star(i,sharedDims); %%% OPTIONAL!!!
        
        % Make the shared dimensions for the current NN the same as the 
        % ones of the closest NN to y_star but from the Z dataset
       % x_cur(sharedDims) = model.X(ind(1),sharedDims); %%%%% OPTIONAL #2 !!!!
        
        %[ZpredMu(k,:), ZpredSigma(k,:)] = vargplvmPosteriorMeanVar(model.comp{infMod}, model.X(ind(k),:));
        YpredMu(k,:) = vargplvmPosteriorMeanVar(sharedModel.comp{1}, x_cur);
        ZpredMu(k,:) = vargplvmPosteriorMeanVar(sharedModel.comp{2}, x_cur);
    end
    fprintf('\n\n');
    
    
    %-- Plots
    
    if exist('makePlots') & ~makePlots
        continue
    end
    % Open a big figure (first 2 args control the position, last 2 control
    % the size)
    figure('Position',[scrsz(3)/100.86 scrsz(4)/6.666 scrsz(3)/2.0457 scrsz(4)/4.5],...
        'Name',['Fig: ' num2str(i) ' (Exp: ' num2str(experimentNo) ')'],'NumberTitle','off')
    numRows = 1;
    numCols = ceil(2*(numberOfNN+1)/numRows);
    plotCounter = 1;

    subplot(numRows, numCols, 1)
    imagesc(reshape(y_star,height,width)), title(['Given']), colormap('gray')
        axis off 

    for k=1:numberOfNN
        subplot(numRows, numCols, k+plotCounter)
        imagesc(reshape(YpredMu(k,:), height, width)), title(['NN: ' num2str(curNNall(k))]), colormap('gray');
        axis off 
    end
    for k=1:numberOfNN
        subplot(numRows, numCols, k+plotCounter+numberOfNN)
        imagesc(reshape(ZpredMu(k,:), height, width)), title(['NN: '  num2str(curNNall(k))]), colormap('gray');
        axis off 
    end
end

