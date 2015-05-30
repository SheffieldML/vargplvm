% DEMYALESVARGPLVM4 Run the Shared Var. GP-LVM on a subset of the Yale
% faces.
% DESC Run the Shared Var. GP-LVM on a subset of the Yale faces. The code
% for creating this subset out of raw images exists in comments. Unlike
% demYaleSvargplvm1, this demo is not a wrapper, it can be used as a
% standalone demo.
%
% VARGPLVM

%{
ca;clear;experimentNo = 1;indTr = 40:200;demSitranFaces1
% TO load:
ca;clear;experimentNo=1;indTr=40:200;trainModel=0;demSitranFaces1

%}


% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);


% Create the dataset out of the images.
baseDir=[localDatasetsDirectorySmall 'sharedVargplvm' filesep 'sitranFaces'];
selDirs = {'andreas','alan','javier'};

for d=1:length(selDirs)
    dirFrom=[baseDir filesep selDirs{d}];
    a=dir(dirFrom);
    counter = 0;
    for i=1:length(a)
        if length(a(i).name)>4
            im = imread([dirFrom filesep a(i).name]);
            im = im(:,:,1);
            %imagesc(im), colormap('gray'); title(a(i).name), pause
            counter = counter+1;
            Yall{d}(counter,:)=im(:)';
        end
    end
    Yall{d} = double(Yall{d});
end
height = size(im,1);
width = size(im,2);

if ~exist('trainModel', 'var'), trainModel = false; end
if ~exist('itNo')         ,  itNo = [500 1000 1500];              end     % Default: 2000
if ~exist('indPoints')    ,  indPoints = 45;          end     % Default: 49
if ~exist('initVardistIters'), initVardistIters = 800;      end
if ~exist('mappingKern')   ,  mappingKern ='rbfardjit'; end

% Set to empty value {} to work with toy data
if ~exist('dataSetNames')    ,    dataSetNames = selDirs;    end
if ~exist('invWidthMult'),       invWidthMult = 5;                     end
if ~exist('dataType'), dataType = 'sitran'; end
if ~exist('latentDimPerModel'), latentDimPerModel = 7; end
if ~exist('experimentNo'), experimentNo = 404; end
if ~exist('doPredictions'), doPredictions = false; end
% If this is true, then the model is in "D > N" mode.
if ~exist('DgtN'), DgtN = true; end
% Create initial X by doing e.g. ppca in the concatenated model.m's or by
% doing ppca in the model.m's separately and concatenate afterwards?
if ~exist('initial_X'), initial_X = 'separately'; end % Other options: 'concatenated'
% Which indices to use for training, rest for test
if ~exist('indTr'), indTr = -1; end

enableParallelism = 0;
lbls = [height width];

if exist('pyramid','var') % extract pyramid representation of the images
    if pyramid
        for e=1:size(Y,2)
            Yall{e} = im2pyramid(Yall{e}, lbls(1), lbls(2), 4);
        end
    end
end

numberOfDatasets = length(Yall);
height = lbls(1); width = lbls(2);


% %{
for i=1:size(Yall{1},1)
    for d=1:numberOfDatasets
        subplot(1,numberOfDatasets,d)
        imagesc(reshape(Yall{d}(i,:),height, width)), colormap('gray');
        axis off
    end
    if i==1
        pause
    else
        pause(0.01)
    end
end
% %}



clear 'd'
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
        X_init = [];
        fprintf('# Initialising X by performing ppca in each observed (scaled) dataset separately and then concatenating...\n');
        for dd = 1:numberOfDatasets
            X_init = [X_init ppcaEmbed(m{dd},latentDimPerModel)];
        end
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
    model = svargplvmRestorePrunedModel(prunedModel, Ytr);
end


%%
axFs = 18;
titleFs = 0;
legFs = 20;

ca
sc = svargplvmShowScales(model,0);
maxScales1 = max(sc{1});   maxScales2 = max(sc{2});   maxScales3 = max(sc{3});
sc{1} = sc{1}./maxScales1; sc{2} = sc{2}./maxScales2; sc{3} = sc{3}./maxScales3;
sc{1}=sigmoid(sc{1}*25)-0.5;
sc{2}=sigmoid(sc{2}*25)-0.5;
sc{3}=sigmoid(sc{3}*25)-0.5;

x=1:size(sc{1},2);
K=1.2;
bar1=bar(x, sc{1}, 'FaceColor', 'b', 'EdgeColor', 'b'); 
set(gca, 'YtickLabel',[]);     set(bar1,'BarWidth',K/1.5);        hold on;
L = get(gca,'XLim');
set(gca,'XTick',1:1:model.q)


bar2=bar(x, sc{2}, 'FaceColor', 'r', 'EdgeColor', 'r'); 
set(gca, 'YtickLabel',[]);     set(bar2,'BarWidth',K/2);

bar3=bar(x, sc{3}, 'FaceColor', 'k', 'EdgeColor', 'k'); 
set(gca, 'YtickLabel',[]);     set(bar3,'BarWidth',K/4);


hold off;
%lg=legend(['scales1'],['scales2']);
%set(lg, 'FontSize',legFs)
set(gca, 'FontSize', axFs);
%%
if ~exist('resultsDynamic')  resultsDynamic = 0; end
if ~resultsDynamic
    return
end

for i=1:length(model.comp)
    if model.comp{i}.DgtN
        model.comp{i}.m = model.comp{i}.mOrig;
    end
end
for i=1:length(model.comp)
    model.comp{i}.vis.startDim = {8,16};
    model.comp{i}.vis.startPos = model.vardist.means(60,:);
end


% ca
% v = 1;
% modelVis = model.comp{v};
% if resultsDynamic
%     %bar(model.comp{v}.kern.comp{1}.inputScales);
%     % The following causes OUTOFMEMORY exception except from when we prune the
%     % video dimensions: (Note: also, lvmVisualise was changed a bit so that no
%     % dynamic slides is presented, because otherwise a strange error occurs).
%     modelVis.y = model.comp{v}.y;
%     %[modelP, newHeight, newWidth] = vargplvmReduceVidModel(modelVis, height, width, 2,2);
%     modelP = modelVis; newHeight = height; newWidth = width;
%     lvmVisualiseGeneral(modelP, [], 'imageMRDVisualise', 'imageMRDModify', false, [newHeight newWidth],0,0,1);
%     figure; svargplvmShowScales(model)
%     clear modelP
% end


%%
if resultsDynamic
    if exist('modelOrig','var'), model = modelOrig; end
    modelOrig = model;
    clear global visualiseInfo
    for v=1:model.numModels
        model.comp{v}.y = Ytr{v};
        reduction = 1; % 4
        opt.showVariance = 1;    opt.showInducing =1 ;
        %opt.showVariance=0; opt.showInducing=0;
        [model.comp{v}, newHeight, newWidth] = vargplvmReduceVidModel(model.comp{v}, height, width, reduction, reduction);
        imv{v}='imageMRDVisualise';
        imm{v}='imageMRDModify';
        args{v}={[newHeight newWidth],0,0,1};
    end
    lvmVisualiseMRD(model,[],imv,imm,opt,args);
    figure;svargplvmShowScales(model);
    set(gca, 'FontSize', 18);
end

%%

%{
%%
ca
pos = 58; 
dim = 17;
nSamp = 15;
ysamp = nan(nSamp, size(modelVis.y,2));
%xsamp = nan(nSamp, size(model.vardist.means,2));
ysamp(1,:) = modelVis.y(pos,:);
%xsamp(1,:) = model.vardist.means(1,:);
minX = min(model.vardist.means(:,dim));
maxX = max(model.vardist.means(:,dim));
rr = abs(maxX - minX);

xx = linspace(3.4,maxX+0.5,nSamp);
xsamp = repmat(model.vardist.means(pos,:), nSamp,1);
xsamp(:,dim) = xx';

for n=1:nSamp
    ysamp(n,:) = vargplvmPosteriorMeanVar(modelVis, xsamp(n,:));
end

pp = 0.05;


%%
ca
imagesc(reshape(ysamp(n,:), height, width));  colormap('gray')
pause
for j=1:10
    for n=2:size(ysamp,1)-1
        imagesc(reshape(ysamp(n,:), height, width));
        colormap('gray')
        pause(pp);
    end
    
    for n=size(ysamp,1):-1:1
        imagesc(reshape(ysamp(n,:), height, width));
        colormap('gray')
        pause(pp);
    end
end
%}

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