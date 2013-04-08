% DEMHUMANPOSETRAINED
% SHEFFIELDML

% E.g: To just see the animation:
% clear;dynUsed=1;showPlots=0;saveFigures=0;testOnTraining=0;numberTestPoints=-1;resultsDynamic=0;showVideo=1;demHumanPoseTrained


if ~exist('showPlots'), showPlots=1; end
if ~exist('saveFigures'), saveFigures=false; end
if ~exist('numberTestPoints'), numberTestPoints = 10; end
% Set to 1 to test on the training data itself, set to 0 to use the test
% dataset.
if ~exist('testOnTraining'), testOnTraining = 1; end
if ~exist('resultsDynamic'),  resultsDynamic = 1; end
if ~exist('showVideo'), showVideo = 0; end

pathToSave = []; %'../../vargplvm/matlab/Results/CVPR/diagrams/humanPose/trainingData/';
addpath(genpath('../'))


if dynUsed
    load demHumanPoseSvargplvm120
    pathToSave = [pathToSave 'dynamic/'];
else
    load demHumanPoseSvargplvm311
    pathToSave = [pathToSave 'static/'];
end

subsample = 2;
seqToKeep = [1 3 4 5 6];
testSeq = 8;

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



%% ----------------------

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


v = infMod;
modelVis = model.comp{v};
if resultsDynamic
    figure
    modelVis.y = Ytr{v};
    modelVis.dynamics = [];
    
   [void, connect] = mocapLoadTextData('run1');   
   dataType = 'xyzankur';
   lvmVisualise(modelVis, [], [dataType 'Visualise'], [dataType 'Modify'],2);
   figure
end
svargplvmShowScales(model);


%%
%---------------------------- PREDICTIONS ---------------


%%%TEMP
%Y_test = Y_test(1:10,:);Z_test = Z_test(1:10,:);
%%%



if dynUsed
    load demHumanPoseSvargplvmPred8Seq120
    model.dynamics.t_star = timeStampsTest;
    model.comp{1}.dynamics.t_star = timeStampsTest;
    model.comp{2}.dynamics.t_star = timeStampsTest; 
    
    for i=1:size(Z_test,1)
        % initialize the latent points using the nearest neighbour from the training data
        dst = dist2(Y_test(i,:), Ytr{1});
        [mind, mini(i)] = min(dst);
    end
    miniNNY = mini;
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
    Xnn = modelOrig.X;
else
    load demHumanPosePred8Seq311
    Xnn = model.X;
    modelOrig = model;
end




if testOnTraining
    perm = randperm(model.N);
    testInd = perm(1:numberTestPoints);
  %  testInd = 1:5:5*numberTestPoints; %%%
  testInd = 1:4:N{1};
else
    Yts{obsMod} = Y_test;
    Yts{infMod} = Z_test;
    perm = randperm(size(Yts{obsMod},1));
    testInd = perm(1:numberTestPoints);
    if numberTestPoints == -1
        testInd = 1:size(Yts{1},1); 
    end
    ZpredAll = zeros(size(Z_test));
    indsAll = zeros(size(Z_test,1),1);   
end

scrsz = get(0,'ScreenSize');

numberOfNN = 4;
for i=1:length(testInd)
    curInd = testInd(i);
    if testOnTraining
        y_star = Ytr{1}(curInd,:);
        z_star = Ytr{2}(curInd,:);
        sil_star = Yim(curInd,:);
        
        x_star(i,:) = model.comp{obsMod}.vardist.means(curInd,:);
        varx_star(i,:) = model.comp{obsMod}.vardist.covars(curInd,:);
        NNstart = 2;
    else
        y_star = Y_test(curInd,:);
        z_star = Z_test(curInd,:);
        sil_star = Yim_test(curInd,:);
        dst = dist2(y_star, model.comp{obsMod}.y);
        [mind, mini] = min(dst);
        [mind,indNN]=sort(dst);
        indNN = indNN(1:numberOfNN);
        miniAll(i) = mini;
        NNstart = 1;
    end
    
    [ind2,distInd] = nn_class(modelOrig.X(:,sharedDims), x_star(i,sharedDims), numberOfNN, 'euclidean');
    indsAllOrig(i) = ind2(1);
        
    nnMode=1; %%%%%%
    if nnMode==1
        w = s1(sharedDims)+s2(sharedDims);
        [ind, distInd] = nn_class2(model.X(:,sharedDims), x_star(i,sharedDims), numberOfNN, 'weighted', w);
    else
        [ind, distInd] = nn_class(Xnn(:,sharedDims), x_star(i,sharedDims), numberOfNN, 'euclidean');
    end
    indsAll(i) = ind(1);
    
    ZpredMu = zeros(length(ind), size(model.comp{infMod}.y,2));
    ZpredSigma = zeros(length(ind), size(model.comp{infMod}.y,2));
    
    
    for k=1:numberOfNN
        x_cur = model.X(ind(k),:);
        %   x_cur(sharedDims) = x_star(i,sharedDims); %%% OPTIONAL!!!
        %---- OPTIONAL 3
          xcurOrig  = x_cur(sharedDims);
          s1new = s1/sum(s1);
          x_cur(sharedDims) = s1new(sharedDims).*x_star(i,sharedDims) + (1-s1new(sharedDims)).*xcurOrig;
        %----
        
        % Optional 4
        % x_cur(sharedDims) = x_star(i, sharedDims);
        % x_cur = svargplvmOptimiseDimVar(model.comp{2}, x_cur, privateDims, 0, 400, 0, 'scg');
        
        %[ZpredMu(k,:), ZpredSigma(k,:)] = vargplvmPosteriorMeanVar(model.comp{infMod}, model.X(ind(k),:));
        if ~testOnTraining
            if dynUsed
                ZpredMu(k,:) = vargplvmPosteriorMeanVar(model.comp{infMod}, x_cur);
            else
                ZpredMu(k,:) = vargplvmPosteriorMeanVar(model.comp{infMod}, x_cur, varx_star(i,:));
            end
        else
            ZpredMu(k,:) = Ytr{2}(ind(k),:);
        end
    end
    ZpredAll(i,:) = ZpredMu(1,:);
    
    if exist('makePlots') & ~makePlots
        continue
    end
    
    %-- Plots
    if showPlots
        
        % Sil: svargplvm
        handle1 = figure;%('Position',[scrsz(3)/100.86 scrsz(4)/6.666 scrsz(3)/3.0457 scrsz(4)/3.0682]);
        numRows = 1;
        if testOnTraining
            numCols = numberOfNN;
        else
            numCols = numberOfNN+1;
        end
        %  ha = tight_subplot(numRows,numCols,[.01 .01],[.01 .01],[.01 .01]);
        ha = tight_subplot(numRows,numCols,[0 0],[0 0],[0 0]);
        axes(ha(1));
        h=imagesc(reshape(sil_star,height,width));colormap('gray'),axis equal; axis tight;axis off; %axis equal; %axis tight; %axis off;
        counter = 2;
        for k=NNstart:numberOfNN
            axes(ha(counter));
            h=imagesc(reshape(Yim(ind(k),:),height, width));colormap('gray'),axis equal; axis tight;axis off;%axis equal; %axis tight; %axis off;
            counter = counter + 1;
        end
        
        % Pose, Svargplvm
        handle2 = figure;
        numRows = 1;
        if testOnTraining
            numCols = numberOfNN;
        else
            numCols = numberOfNN+1;
        end
        %ha = tight_subplot(numRows,numCols,[.01 .01],[.01 .01],[.01 .01]);
        ha = tight_subplot(numRows,numCols,[0 0],[0 0],[0 0]);
        axes(ha(1));
        h=xyzankurVisualise2(z_star);axis off; %axis equal; %axis tight; %axis off;
        counter = 2;
        for k=NNstart:numberOfNN
            axes(ha(counter));
            h=xyzankurVisualise2(ZpredMu(k,:)); axis off; %axis equal; %axis tight; %axis off;
            counter = counter + 1;
        end
        
        %-- To save figures
        if saveFigures
            print( handle1, '-depsc', [pathToSave num2str(curInd) '_sill.eps'])
            print( handle2, '-depsc', [pathToSave num2str(curInd) '_pose.eps'])
            close(handle1)
            close(handle2)
        end
    end
end


%%% !!!!!! SOmething is wrong with the NN results........ % TODO
if dynUsed
    meanPose = repmat(mean(Ytr{2}),size(Z_test,1),1);
    errors.meanPose = xyzankurError(meanPose, Z_test);
    errors.NNYspace = xyzankurError(Ytr{2}(miniNNY,:), Z_test);
    errors.NNXspace = xyzankurError(Ytr{2}(indsAllOrig,:),Z_test);
    errors.svargplvm = xyzankurError(ZpredAll, Z_test);
    fprintf('# Mean Pose Error: %d\n', errors.meanPose)
    fprintf('# NN in the Y space Error: %d\n',errors.NNYspace)
    fprintf('# NN in the X space Error: %d\n',errors.NNXspace)
    fprintf('# Svargplvm Error: %d\n', errors.svargplvm)
    prunedModelUpdated = vargplvmPruneModel(modelUpdated);
   % save(['demHumanPoseSvargplvm' num2str(experimentNo) '.mat'],'barmu','lambda','prunedModelUpdated','errors');
    
   % xyzankurAnimCompareMultipleTEMP(Z_test, {ZpredAll, Ytr{2}(mini,:),
   % Ytr{2}(indsAllOrig,:)},-1,{'Gr. Truth', 'Svargplvm', 'NN_Y','NN_X'});
else
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
        if showVideo
           xyzankurAnimCompareMultipleTEMP(Z_test, {ZpredAll, Ytr{2}(miniAll,:), Ytr{2}(indsAll,:)},-1,{'Gr. Truth', 'Svargplvm', 'NN_Y','NN_X'});
        end
       %lvmVisualise(model.comp{2}, [], 'xyzankurVisualise2', 'xyzankurModify');%%%
    end
end