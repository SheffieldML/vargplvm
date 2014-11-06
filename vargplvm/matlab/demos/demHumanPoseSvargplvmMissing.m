% DEMHUMANPOSESVARGPLVMMISSING Run the shared variational GPLVM on various kinds of data
% allowing missing inputs.
% DESC Run the shared variational GPLVM on various kinds of data
% allowing missing inputs.
%
% COPYRIGHT: Andreas C. Damianou, Carl Henrik Ek,  2011
% SEEALSO : demHumanPoseSvargplvm1
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
if ~exist('vardistCovarsMult'),  vardistCovarsMult=1.3;                  end
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


model = svargplvmRestorePrunedModel(prunedModel, Ytr);
if isfield(model, 'dynamics') & ~isempty(model.dynamics)
    dynUsed=1;
else
    dynUsed = 0;
end
makePlots=0;

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
    startInd = 60;
    missingInd = 1:25;
    Y_testOrig = Y_test;
    Y_test(startInd:end, missingInd) = NaN;

    %model.comp{2}.y(startInd:end, missingInd)=NaN; %%%???????????
    disp('# Initializing latent points...');
    indexPresent = setdiff(1:model.comp{1}.d, missingInd);
        
    Yts{obsMod} = Y_test;
    Yts{infMod} = Z_test;
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
            dst = dist2(y_star(indexPresent), Y_test(:,indexPresent));
                        
            [mind, mini(i)] = min(dst);
            %miniAll(i) = mini;
            Init(i,:) = model.vardist.means(mini(i),:);
            vardistx = vardistCreate(model.comp{obsMod}.vardist.means(mini(i),:), model.q, 'gaussian');
            vardistx.covars = model.comp{obsMod}.vardist.covars(mini(i),:);
            model.comp{obsMod}.vardistx = vardistx;
            display=1;
            iters = 250;
            % Find p(X_* | Y_*) which is approximated by q(X_*)
            [x_star(i,:), varx_star, modelUpdated] = vargplvmOptimisePoint(model.comp{obsMod}, vardistx, y_star, display, iters);%%%
        end
    end
    fprintf('# Predicting images from the NN of X_* ');
 
    for i=1:length(testInd)
        curInd = testInd(i);
        y_star = Yts{obsMod}(curInd,:);
        z_star = Yts{infMod}(curInd,:);
            
        numberOfNN = 1;
        % Now we selected a datapoint X_* by taking into account only the
        % private dimensions for Y. Now, based on the shared dimensions of
        % that, we select the closest (in a NN manner) X from the training
        % data.
      %   w = s1(sharedDims)+s2(sharedDims);
      %   [ind, distInd] = nn_class2(model.X(:,sharedDims), x_star(i,sharedDims), numberOfNN, 'weighted', w);
          [ind, distInd] = nn_class(model.X(:,sharedDims), x_star(i,sharedDims), numberOfNN, 'euclidean');
        indsAll(i) = ind(1);
       
        ZpredMu = zeros(length(ind), size(model.comp{infMod}.y,2));
        ZpredSigma = zeros(length(ind), size(model.comp{infMod}.y,2));
        
        
        % Find p(y_*|x_*) for every x_* found from the NN
        for k=1:numberOfNN
            %fprintf('.');
            x_cur = model.X(ind(k),:);
            % Make the shared dimensions of the current NN the same as the
            % ones of the latent point for y_star
           %  x_cur(sharedDims) = x_star(i,sharedDims); %%% OPTIONAL!!!
            
            % Make the shared dimensions for the current NN the same as the
            % ones of the closest NN to y_star but from the Z dataset
            % x_cur(sharedDims) = model.X(ind(1),sharedDims); %%%%% OPTIONAL #2 !!!!
            
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
        errors.NNYspace = xyzankurError(Ytr{2}(mini,:), Z_test);
        errors.NNXspace = xyzankurError(Ytr{2}(indsAll,:),Z_test);
        errors.svargplvm = xyzankurError(ZpredAll, Z_test);
        fprintf('# Mean Pose Error: %d\n', errors.meanPose)
        fprintf('# NN in the Y space Error: %d\n',errors.NNYspace)
        fprintf('# NN in the X space Error: %d\n',errors.NNXspace)
        fprintf('# Svargplvm Error: %d\n', errors.svargplvm)
        xyzankurAnimCompareMultipleTEMP(Z_test, {ZpredAll, Ytr{2}(mini,:), Ytr{2}(indsAll,:)},-1,{'Gr. Truth', 'Svargplvm', 'NN_Y','NN_X'});
    end
else
    %%% Dynamics
    
    model.dynamics.t_star = timeStampsTest;
    model.comp{1}.dynamics.t_star = timeStampsTest;
    model.comp{2}.dynamics.t_star = timeStampsTest;
    

    for i=1:size(Y_test,1)
        % initialize the latent points using the nearest neighbour from the training data
        dst = dist2(Y_test(i,indexPresent), Ytr{1}(:, indexPresent));
        [mind, mini(i)] = min(dst);
    end
     
    
    indexMissingData = startInd:size(Y_test,1);
    
    vardistx = vardistCreate(model.dynamics.vardist.means(mini,:), model.q, 'gaussian');
    vardistx.covars = 0.2*ones(size(vardistx.covars));
    model.vardistx = vardistx;
    model.comp{1}.vardistx = vardistx;
    model.comp{2}.vardistx = vardistx;
    
    iters = 4000;

    [x, varx] = vargplvmOptimiseSeqDyn(model.comp{obsMod}, vardistx, Y_test, 1, iters);
    % keep the optimized variational parameters
    barmu = x;
    lambda = varx;
    
    % Get the variational means and variacnes for the new test sequcen and
    % update the model to be prepared for prediction
    [x_star, varx_star, modelUpdated] = vargplvmDynamicsUpdateModelTestVar(model.comp{obsMod}, barmu, lambda, Y_test);
   
     
    % Latent variables corresponding to the data  with missing dimensions
    Testmeans = x(indexMissingData, :);
    Testcovars = varx(indexMissingData, :);
    
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
    fprintf('# Predicting images from the NN of X_* ');
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
            x_cur = model.X(ind(k),:);%%%%  % modelOrig.X(ind2(k));
            x_cur(sharedDims) = x_star(i,sharedDims); %%% OPTIONAL!!!
            
                        %---- OPTIONAL 3
          %  xcurOrig  = x_cur(sharedDims);
          %  s1new = s1/sum(s1);
          %  x_cur(sharedDims) = s1new(sharedDims).*x_star(i,sharedDims) + (1-s1new(sharedDims)).*xcurOrig;
            %----  
            
            ZpredMu(k,:) = vargplvmPosteriorMeanVar(model.comp{infMod}, x_cur, varx_star(i,:)); % varx_star needed?
        end
        ZpredAll(i,:) = ZpredMu(1,:);
    end
    fprintf('\n');   
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





