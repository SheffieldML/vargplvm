% DEMHUMANPOSESVARGPLVM1 Run the shared variational GPLVM on the humanEva dataset.
%
% COPYRIGHT: Andreas C. Damianou, Carl Henrik Ek,  2011
% SEEALSO : demHumanPosePrepareData
%
% VARGPLVM


% Human pose data with the whole silhouette
%clear ;close all; experimentNo=404; imageSilhouette=1;initLatent='ppca';
%latentDimPerModel=7;dataSetNames =
%{};toyDataCreate='humanPose';initVardistIters = 380;
%itNo = [500 200 200 200 200 200 200 200 200 200];  indPoints=100;TEMPdemSharedVargplvm


%___

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

% Define constants (in a manner that allows other scripts to parametrize
% this one).
if ~exist('experimentNo') ,  experimentNo = 404;      end
if ~exist('itNo')         ,  itNo = [500 500];              end     % Default: 2000
if ~exist('indPoints')    ,  indPoints = 100;          end     % Default: 49
if ~exist('latentDimPerModel')    ,  latentDimPerModel = 6;          end
% Set to 1 to use dynamics or to 0 to use the standard var-GPLVM
if ~exist('dynUsed')      ,  dynUsed = 0;             end
if ~exist('initVardistIters'), initVardistIters = 180;      end
if ~exist('dynamicKern')   ,     dynamicKern = {'rbf', 'white', 'bias'}; end
if ~exist('mappingKern')   ,  mappingKern = {'rbfard2', 'white', 'bias'}; end


% 0.1 gives around 0.5 init.covars. 1.3 biases towards 0.
if ~exist('vardistCovarsMult'),  vardistCovarsMult=[]; end %1.3; % EDIT from submission.
% Set to empty value {} to work with toy data

if ~exist('invWidthMultDyn'),    invWidthMultDyn = 100;                     end
if ~exist('invWidthMult'),       invWidthMult = 5;                     end
if ~exist('initX'),     initX ='ppca';   end % That's for the dynamics initialisation
if ~exist('dataType'), dataType = 'humanPose'; end
if ~exist('enableParallelism'), enableParallelism = 1; end
if ~exist('DgtN'), DgtN = false; end
% Create initial X by doing e.g. ppca in the concatenated model.m's or by
% doing ppca in the model.m's separately and concatenate afterwards?
if ~exist('initial_X'), initial_X = 'separately'; end % Other options: 'together'
% Which indices to use for training, rest for test
if ~exist('indTr'), indTr = -1; end
if ~exist('usePixels'), usePixels = false; end
if ~exist('parallel', 'var'), parallel = false; end
if ~exist('addBetaPrior','var'), addBetaPrior = false; end
if ~exist('priorScale','var'), priorScale = 1; end
if ~exist('useScaleVal','var'), useScaleVal = true; end

demHumanPosePrepareData


dataSetNames={'silhouette', 'pose'};
% X_init = Xp;
%mappingKern = {'linard2', 'white'};
%mappingKern = {'rbfard2', 'white'};
latentDim = 5; % Anything > 2 and < 10
%xyzankurAnim(Z_test, 3);
numberOfDatasets = length(Yall);



%-- Load datasets
for i=1:numberOfDatasets
    Y = Yall{i};
    
    dims{i} = size(Y,2);
    N{i} = size(Y,1);
    if indTr == -1
        indTr = 1:N{i};
    end
    %t{i} = linspace(0, 2*pi, size(Yall{i}, 1)+1)'; t{i} = t{i}(1:end-1, 1);
    
    indTs = setdiff(1:size(Y,1), indTr);
    Ytr{i} = Y(indTr,:);
    Yts{i} = Y(indTs,:);
    
    d{i} = size(Ytr{i}, 2);
end
% timeStampsTraining = t{1}(indTr,1); %timeStampsTest = t(indTs,1);


t = linspace(0, 2*pi, length(indTr)+1)'; t = t(1:end-1, 1);

% Fix times:
prevSeq = 1;
timeStampsTraining = [];
dt=0.05;
for i=1:length(seq)
    t = ([0:(seq(i)-prevSeq)].*dt)';
    prevSeq = seq(i)+1;
    timeStampsTraining = [timeStampsTraining ;t];
end;

dt=t(2)-t(1);
timeStampsTest = ([0:size(Y_test,1)-1].*dt)';

for i=2:numberOfDatasets
    if N{i} ~= N{i-1}
        error('The number of observations in each dataset must be the same!');
    end
end


%{
for i=1:size(Ytr{1},1)
    %figure
    %subplot(1,2,2)
    %imagesc(reshape(Ytr{1}(i,:), height, width)), colormap('gray');
    %subplot(1,2,1)
    handle = xyzankurVisualise2(Ytr{2}(i,:));
    pause;
    %close
    i
end
xyzankurAnim2(Ytr{2})
%}

%%

%-- Options for the models
for i=1:numberOfDatasets
    % Set up models
    options{i} = vargplvmOptions('dtcvar');
    options{i}.kern = mappingKern; %{'rbfard2', 'bias', 'white'};
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
    if useScaleVal
        options{i}.scaleVal = sqrt(var(Ytr{i}(:)));
    end
    if dynUsed
        optionsDyn{i}.type = 'vargpTime';
        optionsDyn{i}.t=timeStampsTraining;
        optionsDyn{i}.inverseWidth=30;
        optionsDyn{i}.seq = seq;
        % Fill in with default values whatever is not already set
        optionsDyn{i} = vargplvmOptionsDyn(optionsDyn{i});
    end
end

%-------------- INIT LATENT SPACE ---%
for i=1:length(Ytr)
    if usePixels && (i ==1)
        m{1} = Ytr{i};
        medTemp = (max(max(m{1})) + min(min(m{1})))/2;
        ind1 = (m{1} > medTemp);
        ind2 = (m{1} <= medTemp);
        m{1}(ind1) = 1;
        m{1}(ind2) = -1;
    else
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
    %    mAll = [mAll m{i}]; % Concatenation (doesn't work if different sizes)
end




embedFunc = str2func([initX 'Embed']);
if strcmp(initial_X, 'separately')
    fprintf('# Initialising X by performing ppca in each observed (scaled) dataset separately and then concatenating...\n');
    %     X_init{1} = ppcaEmbed(m{1},latentDimPerModel);
    %     X_init{2} = ppcaEmbed(m{2},latentDimPerModel);
    if ~exist('latentDimPerModel') || length(latentDimPerModel)~=2
        latentDimPerModel = [7 3];
    end
    if strcmp(initX, 'vargplvm')
        X_init{1} = embedFunc(m{1},latentDimPerModel(1),[],120,30);
        X_init{2} = embedFunc(m{2},latentDimPerModel(2),[],120,30);
    else
        X_init{1} = embedFunc(m{1},latentDimPerModel(1));
        X_init{2} = embedFunc(m{2},latentDimPerModel(2));
    end
    X_init = [X_init{1} X_init{2}];
    %
    %      [U,V] = pca(m{1},7);
    %      X_init{1} = m{1}*V;
    %      [U,V] = pca(m{2},3);
    %      X_init{2} = m{2}*V;
    %      X_init = [X_init{1} X_init{2}];
    
else
    fprintf('# Initialising X by performing ppca in concatenated observed (scaled) data...\n');
    X_init = embedFunc([m{1} m{2}], latentDimPerModel*2);
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
        optionsDyn{i}.t =timeStampsTraining;
        optionsDyn{i}.inverseWidth=invWidthMultDyn; % Default: 100
        optionsDyn{i}.initX = X_init; % initX; % That's for the dynamics
        
        kern = kernCreate(optionsDyn{i}.t, dynamicKern); % Default: {'rbf','white','bias'}
        
        
        %---- Default values for dynamics kernel
        if strcmp(kern.comp{2}.type, 'white')
            kern.comp{2}.variance = 1e-2; % Usual values: 1e-1, 1e-3
        end
        
        if strcmp(kern.comp{2}.type, 'whitefixed')
            if ~exist('whiteVar')
                whiteVar = 1e-6;
            end
            kern.comp{2}.variance = whiteVar;
            fprintf(1,'# fixedwhite variance: %d\n',whiteVar);
        end
        
        if strcmp(kern.comp{1}.type, 'rbfperiodic')
            if exist('periodicPeriod')
                kern.comp{1}.period = periodicPeriod;
            end
            fprintf(1,'# periodic period: %d\n',kern.comp{1}.period);
        end
        
        % The following is related to the expected number of
        % zero-crossings.(larger inv.width numerator, rougher func)
        if ~strcmp(kern.comp{1}.type,'ou')
            kern.comp{1}.inverseWidth = optionsDyn{i}.inverseWidth./(((max(t)-min(t))).^2);
            kern.comp{1}.variance = 1;
        end
        %------
        
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



%modelInit = model;%%%TEMP

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


%-- Define what level of parallelism to use (w.r.t submodels or/and w.r.t
% datapoints).
%{
fprintf('# Parallel computations w.r.t the submodels!\n');
model.parallel = 1;
model = svargplvmPropagateField(model,'parallel', 1);
%}
if parallel
    fprintf('# Parallel computations w.r.t the datapoints!\n');
    model.vardist.parallel = 1;
    for i=1:model.numModels
        model.comp{i}.vardist.parallel = 1;
    end
end
% %}

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


if addBetaPrior
    %--- Config: where and what prior to add
    meanSNR = 150;                       % Where I want the expected value of my inv gamma if it was on SNR
    priorName = 'invgamma';              % What type of prior
    scale = priorScale*model.N;    % 'strength' of prior.
    %----
    for i=1:length(model.comp)
        if isfield(model.comp{i}, 'mOrig'), varData = var(model.comp{i}.mOrig(:));  else,  varData = var(model.comp{i}.m(:));  end
        meanB = meanSNR./varData;
        a=0.08;%1.0001; % Relatively large right-tail
        b=meanB*(a+1); % Because mode = b/(a-1)
        % Add the prior on parameter 'beta'. The prior is specified by its name
        % (priorName) and its parameters ([a,b])
        model.comp{i} = vargplvmAddParamPrior(model.comp{i}, 'beta', priorName, [a b]);
        
        % Add a scale to the prior ("strength") and save this version of the model.
        model.comp{i}.paramPriors{1}.prior.scale = scale;
    end
    params = modelExtractParam(model);
    model = modelExpandParam(model, params);
end

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

%%

if ~exist('resultsDynamic')  resultsDynamic = 0; end

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
    %lvmVisualise(modelP, [], 'imageVisualise', 'imageModify', [newHeight newWidth],0,0,1);
    lvmVisualiseGeneral(modelP, [], 'imageMRDVisualise', 'imageMRDModify', false, [newHeight newWidth],0,0,1);
    clear modelP
end




if ~exist('doPredictions'), doPredictions = 0; end
if ~doPredictions
    return
end

%%
%---------------------------- PREDICTIONS ---------------

% Set to 1 to test on the training data itself, set to 0 to use the test
% dataset.
if ~exist('testOnTraining')
    testOnTraining=0;
end

% 1 is for the HoG image features. 2 is for the pose features.
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
numberTestPoints = 10;
if testOnTraining
    perm = randperm(model.N);
    testInd = perm(1:numberTestPoints);
else
    Yts{obsMod} = Y_test;
    Yts{infMod} = Z_test;
    perm = randperm(size(Yts{obsMod},1));
   % testInd = perm(1:numberTestPoints);
    testInd = 1:size(Yts{1},1);
    ZpredAll = zeros(size(Z_test));
    indsAll = zeros(size(Z_test,1),1);
end

scrsz = get(0,'ScreenSize');

x_star = zeros(length(testInd), size(model.X,2));
if ~(~testOnTraining & dynUsed) % If we have a test set and dynamics, then we need a different inference procedure
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
            miniAll(i) = mini;
            Init(i,:) = model.vardist.means(mini,:);
            vardistx = vardistCreate(model.comp{obsMod}.vardist.means(mini,:), model.q, 'gaussian');
            vardistx.covars = model.comp{obsMod}.vardist.covars(mini,:);
            model.comp{obsMod}.vardistx = vardistx;
            display=1;
            iters = 250;
            % Find p(X_* | Y_*) which is approximated by q(X_*)
            [x_star(i,:), varx_star(i,:), modelUpdated] = vargplvmOptimisePoint(model.comp{obsMod}, vardistx, y_star, display, iters);%%%
        end
    end
    
    fprintf('# Predicting images from the NN of X_* ');
    for i=1:length(testInd)
        curInd = testInd(i);
        numberOfNN = 9;
        if ~testOnTraining
            y_star = Yts{obsMod}(curInd,:);
            z_star = Yts{infMod}(curInd,:);
            numberOfNN=1;
        end
        % Now we selected a datapoint X_* by taking into account only the
        % private dimensions for Y. Now, based on the shared dimensions of
        % that, we select the closest (in a NN manner) X from the training data.
        % fprintf('# Finding the %d NN of X_* with the training X based only on the shared dims.\n', numberOfNN);
        % w = s1(sharedDims)+s2(sharedDims);
        % [ind, distInd] = nn_class2(model.X(:,sharedDims), x_star(i,sharedDims), numberOfNN, 'weighted', w);
        [ind, distInd] = nn_class(model.X(:,sharedDims), x_star(i,sharedDims), numberOfNN, 'euclidean');
        % [ind, distInd] = nn_class(model.X, x_star(i,:), numberOfNN, 'euclidean');%%%%?
        indsAll(i) = ind(1);
        
        ZpredMu = zeros(length(ind), size(model.comp{infMod}.y,2));
        ZpredSigma = zeros(length(ind), size(model.comp{infMod}.y,2));
        
        
        % Find p(y_*|x_*) for every x_* found from the NN
        for k=1:numberOfNN
            x_cur = model.X(ind(k),:);
            % Make the shared dimensions of the current NN the same as the
            % ones of the latent point for y_star
            %   x_cur(sharedDims) = x_star(i,sharedDims); %%% OPTIONAL!!!
            
            % Make the shared dimensions for the current NN the same as the
            % ones of the closest NN to y_star but from the Z dataset
            % x_cur(sharedDims) = model.X(ind(1),sharedDims); %%%%% OPTIONAL #2 !!!!
            
            %---- OPTIONAL 3
            xcurOrig  = x_cur(sharedDims);
            s1new = s1/sum(s1);
            x_cur(sharedDims) = s1new(sharedDims).*x_star(i,sharedDims) + (1-s1new(sharedDims)).*xcurOrig;
            %----
            
            % Optional 4
            % x_cur(sharedDims) = x_star(i, sharedDims);
            % x_cur = svargplvmOptimiseDimVar(model.comp{2}, x_cur, privateDims, 0, 400, 0, 'scg');
            
            
            %[ZpredMu(k,:), ZpredSigma(k,:)] = vargplvmPosteriorMeanVar(model.comp{infMod}, model.X(ind(k),:));
            %if ~testOnTraining
            ZpredMu(k,:) = vargplvmPosteriorMeanVar(model.comp{infMod}, x_cur);
            %else
            %    ZpredMu(k,:) = model.comp{infMod}(ind(k),:);
            %end
        end
        ZpredAll(i,:) = ZpredMu(1,:);
        
        %-- Plots
        if exist('makePlots') & ~makePlots
            continue
        end
        % Open a big figure (first 2 args control the position, last 2 control
        % the size)
        figure('Position',[scrsz(3)/100.86 scrsz(4)/6.666 scrsz(3)/1.0457 scrsz(4)/1.0682],...
            'Name',['Fig: ' num2str(i) ' (Exp: ' num2str(experimentNo) ')'],'NumberTitle','off')
        numRows = 3;
        
        if testOnTraining
            numCols = ceil((numberOfNN+1)/numRows)*2;
            plotCounter = 2;
        else
            % For the real test image!
            numCols = ceil((numberOfNN+2)/numRows);
            plotCounter = 2;
        end
        
        
        if ~testOnTraining
            sil_star = Yim_test(curInd,:);
            pose_star = Z_test(curInd,:);
            subplot(numRows, numCols, 1)
            imagesc(reshape(sil_star,height,width))
            if obsMod == 1
                title(['Given (image #' num2str(curInd) ')']), colormap('gray')
            else
                title('Corresponding');
            end
            subplot(numRows, numCols, 2)
            handle = xyzankurVisualise2(pose_star);
            if obsMod == 1
                title(['Given (image #' num2str(curInd) ')']), colormap('gray')
            else
                title('Corresponding');
            end
            
            for k=1:numberOfNN
                subplot(numRows, numCols, k+plotCounter)
                if infMod == 1 &&  exist('imageSilhouette') && imageSilhouette
                    imagesc(reshape( ZpredMu(k,:),height,width)), title(['NN #' num2str(k)]), colormap('gray')
                else
                    handle = xyzankurVisualise2(ZpredMu(k,:)); title(['NN #' num2str(k)])
                end
            end
            %mserror(i) = mean(abs(ZpredMu(1,:) - Z_test(i,:)));
        else
            subplot(numRows, numCols, 1)
            imagesc(reshape(Yim(curInd,:),height,width)), title(['Original y (image #' num2str(curInd) ')']), colormap('gray')
            subplot(numRows, numCols, 2)
            handle = xyzankurVisualise2(model.comp{2}.y(curInd,:));
            % Start plotting from 2, the first is always the same as
            for k=2:numberOfNN
                % Start from k=2, the first NN we know it's the same as the
                % given point.
                subplot(numRows, numCols, k+plotCounter-1)
                imagesc(reshape(Yim(ind(k),:),height,width)), title(['NN #' num2str(k)]), colormap('gray')
                subplot(numRows, numCols, k+plotCounter)
                handle = xyzankurVisualise2(model.comp{2}.y(ind(k),:)); title(['NN #' num2str(k)])
                plotCounter = plotCounter+1;
            end
        end
        %     pause
        %     close
    end
    if ~testOnTraining
        meanPose = repmat(mean(Ytr{2}),size(Y_test,1),1);
        errors.meanPose = xyzankurError(meanPose, Z_test);
        errors.NNYspace = xyzankurError(Ytr{2}(miniAll,:), Z_test);
        errors.NNXspace = xyzankurError(Ytr{2}(indsAll,:),Z_test);
        errors.svargplvm = xyzankurError(ZpredAll, Z_test);
        fprintf('# Mean Pose Error: %d\n', errors.meanPose)
        fprintf('# NN in the Y space Error: %d\n',errors.NNYspace)
        fprintf('# NN in the X space Error: %d\n',errors.NNXspace)
        fprintf('# Svargplvm Error: %d\n', errors.svargplvm)
        % xyzankurAnimCompareMultipleTEMP(Z_test, {ZpredAll, Ytr{2}(miniAll,:), Ytr{2}(indsAll,:)},-1,{'Gr. Truth', 'Svargplvm', 'NN_Y','NN_X'});
    end
else
    %%% Dynamics
    
    model.dynamics.t_star = timeStampsTest;
    model.comp{1}.dynamics.t_star = timeStampsTest;
    model.comp{2}.dynamics.t_star = timeStampsTest;
    
    for i=1:size(Z_test,1)
        % initialize the latent points using the nearest neighbour from the training data
        dst = dist2(Y_test(i,:), Ytr{1});
        [mind, mini(i)] = min(dst);
    end
    
    vardistx = vardistCreate(model.dynamics.vardist.means(mini,:), model.q, 'gaussian');
    vardistx.covars = 0.2*ones(size(vardistx.covars));
    model.vardistx = vardistx;
    model.comp{1}.vardistx = vardistx;
    model.comp{2}.vardistx = vardistx;
    iters = 4500;
    % Do also reconstruction in test data
    [x, varx] = vargplvmOptimiseSeqDyn(model.comp{obsMod}, vardistx, Y_test, 1, iters);
    % keep the optimized variational parameters
    barmu = x;
    lambda = varx;
    
    % Get the variational means and variacnes for the new test sequcen and
    % update the model to be prepared for prediction
    [x_star, varx_star, modelUpdated] = vargplvmDynamicsUpdateModelTestVar(model.comp{obsMod}, barmu, lambda, Y_test);
    modelOrig = model;
    model.comp{obsMod} = modelUpdated;
    model.vardist = modelUpdated.vardist;
    model.comp{1}.vardist = modelUpdated.vardist;
    model.comp{2}.vardist = modelUpdated.vardist;
    model.dynamics = modelUpdated.dynamics;
    model.comp{1}.dynamics = modelUpdated.dynamics;
    model.comp{2}.dynamics = modelUpdated.dynamics;
    model.X = model.vardist.means;
    model.comp{1}.X = model.X;
    model.comp{2}.X = model.X;
    
    numberOfNN = 2;
    
    fprintf('# Finding the %d NN of X_* with the training X based only on the shared dims.\n', numberOfNN);
    ZpredAll = zeros(size(Z_test));
    indsAll = zeros(size(Z_test,1),1);
    indsAllOrig = zeros(size(Z_test,1),1);
    
    for i=1:size(Z_test,1)
        [ind2,distInd] = nn_class(modelOrig.X(:,sharedDims), x_star(i,sharedDims), numberOfNN, 'euclidean');
        indsAllOrig(i) = ind2(1);
        
        % Actually this is a weighted NN.
        %  w = s1(sharedDims)+s2(sharedDims);
        %  [ind, distInd] = nn_class2(model.X(:,sharedDims), x_star(i,sharedDims), numberOfNN, 'weighted', w);
        [ind, distInd] = nn_class(model.X(:,sharedDims), x_star(i,sharedDims), numberOfNN, 'euclidean');
        indsAll(i) = ind(1);
        ZpredMu = zeros(length(ind), size(model.comp{infMod}.y,2));
        ZpredSigma = zeros(length(ind), size(model.comp{infMod}.y,2));
        % Find p(y_*|x_*) for every x_* found from the NN
        for k=1:numberOfNN
            % fprintf('.');
            x_cur = model.X(ind(k),:);%%%%  % modelOrig.X(ind2(k));
            x_cur(sharedDims) = x_star(sharedDims); %%% OPTIONAL!!!
            
            %---- OPTIONAL 3
            %  xcurOrig  = x_cur(sharedDims);
            %  s1new = s1/sum(s1);
            %  x_cur(sharedDims) = s1new(sharedDims).*x_star(i,sharedDims) + (1-s1new(sharedDims)).*xcurOrig;
            %----
            
            % Optional 4
            % x_cur(sharedDims) = x_star(i, sharedDims);
            % x_cur = svargplvmOptimiseDimVar(model.comp{2}, x_cur, privateDims, 0, 400, 0, 'scg');
            
            
            ZpredMu(k,:) = vargplvmPosteriorMeanVar(model.comp{infMod}, x_cur);%, varx_star(i,:)); % varx_star needed?
        end
        ZpredAll(i,:) = ZpredMu(1,:);
    end
    fprintf('\n\n');
    meanPose = repmat(mean(Ytr{2}),size(Z_test,1),1);
    errors.meanPose = xyzankurError(meanPose, Z_test);
    errors.NNYspace = xyzankurError(Ytr{2}(mini,:), Z_test);
    errors.NNXspace = xyzankurError(Ytr{2}(indsAllOrig,:),Z_test);
    errors.svargplvm = xyzankurError(ZpredAll, Z_test);
    fprintf('# Mean Pose Error: %d\n', errors.meanPose)
    fprintf('# NN in the Y space Error: %d\n',errors.NNYspace)
    fprintf('# NN in the X space Error: %d\n',errors.NNXspace)
    fprintf('# Svargplvm Error: %d\n', errors.svargplvm)
    prunedModelUpdated = vargplvmPruneModel(modelUpdated);
    % save(['demHumanPoseSvargplvm' num2str(experimentNo) '.mat'],'barmu','lambda','prunedModelUpdated','errors');
    
    % xyzankurAnimCompareMultipleTEMP(Z_test, {ZpredAll, Ytr{2}(mini,:), Ytr{2}(indsAllOrig,:)},-1,{'Gr. Truth', 'Svargplvm', 'NN_Y','NN_X'});
end

% plotsDynamic; %%% UNCOMMENT to produce plots

% NN in the pose space
for i=1:158
    dst = dist2(Z_test(i,:), Ytr{2});
    [mind, miniTEMP] = min(dst);
    miniAll2(i) = miniTEMP;
end
xyzankurError(Z_test, Ytr{2}(miniAll2,:))


xyzankurAnimCompareMultipleTEMP(Z_test, {Ytr{2}(mini,:), Ytr{2}(miniAll2,:)}, -1, {'True', 'NN sil space', 'NN pose space'})
