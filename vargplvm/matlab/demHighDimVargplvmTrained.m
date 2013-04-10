% DEMHIGHDIMVARGPLVMTRAINED Perform tasks (predictions, sampling) for an already trained model on high dimensinoal video datasets.
% DESC: Load an already trained model on a video dataset and do sampling or predictions with that.
% !!! The parameters that concern the dataset and the experimentNo must be exactly the same as those for the trained model.
%
% COPYRIGHT: Andreas C. Damianou, Michalis K. Titsias, 2010 - 2011
% SEEALSO: demHighDimVargplvm3.m
% VARGPLVM

%---- This file is in essence the same as demHighDimVargplvm3.m with the training phase omitted. ---
%-------------
% VARGPLVM



% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);


%-- Arguments needed in any case:
%trainedModel, dataSetName, experimentNo, dataSetSplit, <blockSize>, <reconstrIters>

%-- Example 1: white var. 0.06, matern, blocks of points
% trainedModel = 'demOceanVargplvm9';
% dataSetName = 'ocean';
% experimentNo=22;
% dataSetSplit = 'blocks';
% blockSize = 5;
% reconstrIters = 1000; % if predicting with missing inputs

%- optional:
% %testOnTrainingData=1; % optional
% %testOnTrainingTimes=1; % optional
% %dataToKeep %optional
% %noiselessSamples=1; % optional
% %futurePred=times_ahead; % optional
% %superSampling=superSamplingRate; % optional: dt/superSamplingRate
%-


%-- Example 2:  white var.8e-6 , periodic kernel, consecutive points
% trainedModel = 'demOceanVargplvm11';
% dataSetName = 'ocean';
% experimentNo=11;
% dataSetSplit='halfAndHalf';
%-


if ~exist('predWithMs') ,  predWithMs = 0;      end
if ~exist('timesPrediction') ,  timesPrediction = 1;      end
if ~exist('doSampling') ,  doSampling = 0;      end
if ~exist('noiselessSamples'), noiselessSamples=1; end


if ~exist('trainedModel')
    capName = dataSetName;
    capName(1) = upper(capName(1));
    modelType = 'Vargplvm';
    trainedModel = ['dem' capName modelType num2str(experimentNo)];
end

fprintf('# Loading %s: \n', trainedModel);
load(trainedModel);


indPoints=prunedModel.k;
latentDim=prunedModel.q;
for i=1:numel(prunedModel.dynamics.kern.comp)
    dynamicKern{i} = prunedModel.dynamics.kern.comp{i}.type;
end
dynUsed = 1;


% Define constants (in a manner that allows other scripts to parametrize
% this one).
if ~exist('experimentNo') ,  experimentNo = 404;      end
% Set to 1 to keep only the dimensions modelling the head (missa dataset)
if ~exist('cropVideo')    ,  cropVideo = 0;           end
% Set to 1 to remove the frames with the translation (missa dataset)
if ~exist('removeTransl') ,  removeTransl = 0;        end
% Set to 1 to tie the inducing points with the latent vars. X
if ~exist('fixInd')        ,     fixInd = 0;                             end
if ~exist('reconstrIters') ,     reconstrIters = 1000;                   end
% 0.1 gives around 0.5 init.covars. 1.3 biases towards 0.
if ~exist('dataSetSplit')   ,    dataSetSplit = 'everyTwo';              end
if ~exist('blockSize')      ,    blockSize = 8;                          end
% 720x1280
if ~exist('dataSetName')    ,    dataSetName = 'ocean';                  end
if ~exist('testReoptimise') ,    testReoptimise = 1;                     end

if strcmp(dataSetName, 'missa') & strcmp(dataSetSplit,'randomBlocks')
    rand; % TMP (just to make the random seed more convenient! (this results in Ytr close to Yts).
end
if strcmp(dataSetName, 'ocean') & strcmp(dataSetSplit,'randomBlocks')
    for i=1:5, rand; end % TMP (just to make the random seed more convenient! (this results in Ytr close to Yts).
end

% dataSetName = 'susie';
% load susie_352_240

fprintf(1,'\n#----------------------------------------------------\n');
fprintf(1,'# Dataset: %s\n',dataSetName);
fprintf(1,'# ExperimentNo: %d\n', experimentNo);
fprintf(1,'# Latent dimensions: %d\n',latentDim);
fprintf(1,'# Reconstruction iterations: %d\n',reconstrIters);
fprintf(1,'# Tie Inducing points: %d\n',fixInd);
fprintf(1,'# Reoptimise inducing points (for reconstr): %d \n',testReoptimise);
fprintf(1,'# Dynamics used: %d\n', dynUsed);
if dynUsed
    fprintf(1,'# Dynamics kern: ');
    disp(dynamicKern);
end

fprintf(1,'# CropVideo / removeTranslation: %d / %d \n', cropVideo, removeTransl);
if exist('dataToKeep')  fprintf(1,'# DataToKeep: %d \n',dataToKeep); end
%fprintf(1,'#----------------------------------------------------\n');

switch dataSetName
    case 'missa'
        try
            load 'miss-americaHD'
        catch
            Y = vargplvmLoadData(dataSetName);
        end
        width=360;
        height=288;
    case 'ocean'
        try
            load 'DATA_Ocean'
        catch
            Y=vargplvmLoadData('ocean');
        end
        width=1280;
        height=720;
    case 'horse'
        try
            load 'DATA_Horse'
        catch
            Y=vargplvmLoadData('horse');
        end
        %%%%
        Y = Y(90:end,:);
        %%%%
        width=249;
        height=187;
    case 'horse2'
        Y=vargplvmLoadData('horse2');
        width=320;
        height=240;
    case 'horse2cropped'
        Y=vargplvmLoadData('horse2cropped');
        width=249;
        height=187;
    otherwise
        try
            [Y, lbls] = vargplvmLoadData(dataSetName);
            height = lbls(1);
            width = lbls(2);
        catch
            load(dataSetName);
        end
end


if exist('dataToKeep')
    Y = Y(1:dataToKeep,:);
end

% Remove the translation part (only for the "missa" dataset)
if removeTransl
    % For this dataset there is a translation in space between frames 66-103
    if strcmp(dataSetName,'missa')
        Y = [Y(1:64,:) ; Y(103:end,:)];
    else
        error('removeTranslation option only for the missa dataset');
    end
end

% Crop video
if cropVideo
    switch dataSetName
        case 'missa'
            % Keep only the head
            newWidth=144;
            newHeight=122;
            Ynew = zeros(size(Y,1), newWidth*newHeight);
            for i=1:size(Y,1)
                fr=reshape(Y(i,:),height,width);
                fr=imcrop(fr,[129.5 91.5 newWidth-1 newHeight-1]);
                Ynew(i,:) = reshape(fr,1,newWidth*newHeight);
            end
            Y = Ynew; width=newWidth; height=newHeight;
        case 'ocean'
            % Keep only the water
            %[1.5 288.5 1279 432]
            newWidth=width;
            newHeight=round(height/1.6666);
            Ynew = zeros(size(Y,1), newWidth*newHeight);
            for i=1:size(Y,1)
                fr=reshape(Y(i,:),height,width);
                fr=imcrop(fr,[1 288.5 newWidth-1 newHeight-1]);
                Ynew(i,:) = reshape(fr,1,newWidth*newHeight);
            end
            Y = Ynew; width=newWidth; height=newHeight;
        case 'horse'
            % Step 1
            cropVector=[24.5 1.5 224 height]; % horse
            %%%
            newHeight=cropVector(4)-1;
            newWidth=cropVector(3)+1;
            Ynew = zeros(size(Y,1), newWidth*newHeight);
            for i=1:size(Y,1)
                fr=reshape(Y(i,:),height,width);
                fr=imcrop(fr,cropVector);
                Ynew(i,:) = reshape(fr,1,newWidth*newHeight);
            end
            Y = Ynew; width=newWidth; height=newHeight;
            
            % Step 2
            cropVector = [0 0 188 159];
            newHeight=cropVector(4);
            newWidth=cropVector(3);
            Ynew = zeros(size(Y,1), newWidth*newHeight);
            for i=1:size(Y,1)
                fr=reshape(Y(i,:),height,width);
                fr=imcrop(fr,cropVector);
                Ynew(i,:) = reshape(fr,1,newWidth*newHeight);
            end
            Y = Ynew; width=newWidth; height=newHeight;
        case 'dog2'
            cropVector=[11.5 3.5 316 357];
            newHeight=cropVector(4);
            newWidth=cropVector(3)+1;
            Ynew = zeros(size(Y,1), newWidth*newHeight);
            for i=1:size(Y,1)
                fr=reshape(Y(i,:),height,width);
                fr=imcrop(fr,cropVector);
                Ynew(i,:) = reshape(fr,1,newWidth*newHeight);
            end
            Y = Ynew; width=newWidth; height=newHeight;
    end
end

dims = size(Y,2);

trainModel = 1;% Set it to 1 to retrain the model. Set it to 0 to load an already trained one.

% Take a downsample version of the video for training and test on the
% remaining frames
N = size(Y,1);

fprintf(1,'# Preparing the dataset...\n');
switch dataSetSplit
    case 'everyTwo'
        % Training: 1,3,5... Test: 2,4,6...
        % assumes the number of data is even
        Nstar = N/2; Nts = N/2; Ntr = N/2;
        indTr = 1:2:N; indTs = 2:2:N;
    case 'halfAndHalf'
        % Training: 1,2,...N/2  Test: N/2+1, ...,N
        Nstar = round(N/2);
        indTr = 1:size(Y,1)-Nstar; indTs = size(Y,1)-Nstar+1:size(Y,1);
    case 'blocks'
        % Training: 1,2,..B, 2*B+1,.. 3*B,...
        % Test:B+1,...2*B,...   i.e. like "everyTwo" but with blocks
        lastBlockSize = mod(N,blockSize*2);
        mask = [ones(1,blockSize) zeros(1,blockSize)];
        mask = repmat(mask, 1, floor(floor(N/blockSize)/2));
        %md=mod(lastBlockSize,blockSize);
        %md2=lastBlockSize-md;
        %if mask(end)
        %    mask = [mask zeros(1,md2) ones(1,md)];
        %else
        %    mask = [mask ones(1, md2) zeros(1,md)];
        %end
        mask = [mask ones(1,lastBlockSize)]; % always end with the tr. set
        indTr = find(mask);
        indTs = find(~mask);
        if exist('msk')
            indTr = sort([indTr msk]);
            indTs=setdiff(1:N,indTr);
        end
        Nstar = size(indTs,2);
    case 'randomBlocks'
        mask = [];
        lastTrPts = 5;
        r=1; % start with tr. set
        while length(mask)<size(Y,1)-lastTrPts %The last lastTrPts will be from YTr necessarily
            blockSize = randperm(8);
            blockSize = blockSize(1);
            pts = min(blockSize, size(Y,1)-lastTrPts - length(mask));
            if r
                mask = [mask ones(1,pts)];
            else
                mask = [mask zeros(1,pts)];
            end
            r = ~r; % alternate between tr. and test set
        end
        mask = [mask ones(1,lastTrPts)];
        indTr = find(mask);
        indTs = find(~mask);
        if sum(sort([indTr indTs]) - (1:size(Y,1)))
            error('Something went wrong in the dataset splitting...');
        end
        % indPoints = length(indTr); %%%%% temp
        Nstar = length(indTs);
end

if indPoints == -1
    indPoints = length(indTr);
end

Ytr = Y(indTr,:); Yts = Y(indTs,:);
t = linspace(0, 2*pi, size(Y, 1)+1)'; t = t(1:end-1, 1);
timeStampsTraining = t(indTr,1); timeStampsTest = t(indTs,1);


%-- For DEBUG
if exist('testOnTrainingTimes') && testOnTrainingTimes
    timeStampsTest = timeStampsTraining;
end
if exist('testOnTrainingData') && testOnTrainingData
    Yts = Ytr(1:length(timeStampsTest),:);
end
if exist('testOnReverseTrData') && testOnReverseTrData
    Yts = Ytr(end:-1:1,:); timeStampsTest = timeStampsTest(1:size(Yts,1));
end
%--

YtsOriginal = Yts;

fprintf(1,'# Inducing points: %d\n',indPoints);
fprintf(1,'# Dataset size used (train/test) : %d / %d \n', size(Ytr,1), size(Yts,1));
fprintf(1,'# Dataset Split: %s ',dataSetSplit);
if strcmp(dataSetSplit,'blocks'),    fprintf(1,' (blockSize:%d)',blockSize); end
fprintf(1,'\n');
%fprintf(1,'# CropVideo / removeTranslation: %d / %d \n', cropVideo, removeTransl);
fprintf(1,'#----------------------------------------------------\n');

clear Y % Free up some memory

%%
fprintf('# Restoring pruned models...\n');
model = vargplvmRestorePrunedModel(prunedModelTr, Ytr);
%model = vargplvmRestorePrunedModel(prunedModel, Ytr);
%modelInit = vargplvmRestorePrunedModel(prunedModelInit, Ytr); modelTr = vargplvmRestorePrunedModel(prunedModelTr, Ytr);
if testReoptimise
    model.dynamics.reoptimise=1; modelInit.dynamics.reoptimise = 1; modelTr.dynamics.reoptimise = 1;
else
    model.dynamics.reoptimise=0; modelInit.dynamics.reoptimise=0; modelTr.dynamics.reoptimise=0;
end
capName = dataSetName;
capName(1) = upper(capName(1));
modelType = model.type;
modelType(1) = upper(modelType(1));
fileToSave = ['dem' capName modelType num2str(experimentNo) '.mat'];
modelTr.dynamics.t_star = timeStampsTest; model.dynamics.t_star = timeStampsTest;
%%




%%%-----------------------   RECONSTRUCTION --------------------%%%

%model = modelTr;
%clear modelTr %%
model.y = Ytr;
model.m= gpComputeM(model); %%%
model.y=[];
Nstar = size(YtsOriginal,1);
%%% NOTE: If you are reloading and the restoring the "modelUpdated" you will need to:
% modelUpdated = vargplvmRestorePrunedModel(prunedModelUpdated,Ytr,1);
%  modelUpdated.m = modelUpdated.mOrig;
%  modelUpdated.P = modelUpdated.P1 * (modelUpdated.Psi1' * modelUpdated.m);
%%

if exist('futurePred') % futurePred is the number of frames to perdict in future
    dt=timeStampsTraining(end) - timeStampsTraining(end-1);
    if exist('superSampling')
        dt = dt./superSampling;
    end
    if length(futurePred) == 1 % one number, assuming start point the last training
        t_star = timeStampsTraining(end)+dt:dt:timeStampsTraining(end)+dt*futurePred;
    else
        t_star = futurePred; % give the whole training times interval
    end
    %t_star = timeStampsTraining(end):dt:timeStampsTraining(end)+dt*futurePred;
    model.dynamics.t_star = t_star';
    modelTr.dynamics.t_star = t_star';
    timeStampsTest = t_star';
    Nstar = length(timeStampsTest);
    YtsOriginal = zeros(Nstar, model.d);
end

if timesPrediction
    fprintf('# Only times prediction...\n');
    % Prediction using the only information in the test time points
    [Testmeans2 Testcovars2] = vargplvmPredictPoint(model.dynamics, timeStampsTest);
    Varmu2 = vargplvmPosteriorMeanVar(model, Testmeans2, Testcovars2);
    
    if ~exist('futurePred')
        % Mean absolute error per pixel
        errorOnlyTimes = mean(abs(Varmu2(:) - YtsOriginal(:)));
        fprintf(1,'# Error GPLVM: %d\n', errorOnlyTimes);
    end
end


%%


% if you like you can also do prediction after training with partially
% observed test images (also good for de-bugging because this should be better
% than the onlyTimes prediction)



if predWithMs
    %
    display = 1;
    indexP = [];
    Init = [];
    
    % w=360; % width
    % h=288; % height
    w=width; h=height;
    
    
    if ~exist('cut')
        cut = 'vertically';
        if strcmp(dataSetName,'missa') || strcmp(dataSetName,'ADN') || strcmp(dataSetName,'claire')
            cut = 'horizontally';
        end
    end
    
    switch cut
        case 'horizontally'
            if strcmp(dataSetName,'ADN')
                cutPoint=round(h/1.55);
            elseif strcmp(dataSetName,'claire')
                cutPoint=round(h/2);
            elseif strcmp(dataSetName, 'grandma')
                cutPoint=round(h/1.68);
            else %missa
                 cutPoint=round(h/2)+13;
            end
            if exist('cutP')    cutPoint = cutP; end
            mask = [ones(1,cutPoint) zeros(1,h-cutPoint)];
            mask=repmat(mask, 1,w);
            indexMissing = find(mask);
        case 'vertically'
            if strcmp(dataSetName,'missa')
                indexMissing=1:round(size(Yts,2)/1.8);
            elseif strcmp(dataSetName,'dog')
                indexMissing = 1:round(size(Yts,2)/1.70);
            elseif strcmp(dataSetName,'ADN')
                indexMissing = 1:round(size(Yts,2)/1.7);
           % elseif strcmp(dataSetName,'ocean')
           %     indexMissing = 1:round(size(Yts,2)/1.6);
            else
                indexMissing=1:round(size(Yts,2)/2);
            end
    end
    indexPresent = setdiff(1:model.d, indexMissing);
    
    %
    Yts = YtsOriginal;
    
    mini =[];
    for i=1:size(Yts,1)
        % initialize the latent points using the nearest neighbour
        % from he training data
        dst = dist2(Yts(i,indexPresent), Ytr(:,indexPresent));
        % mind: smaller distance %mini(i): smaller index
        [mind, mini(i)] = min(dst);
    end
    
    
    vardistx = vardistCreate(model.dynamics.vardist.means(mini,:), model.q, 'gaussian');
    vardistx.covars = model.dynamics.vardist.covars(mini,:);%0.2*ones(size(vardistx.covars));
    model.vardistx = vardistx;
    
    Yts(:,indexMissing) = NaN;
    model.dynamics.t_star = timeStampsTest;
    iters=reconstrIters;
    
    fprintf(1, '# Partial reconstruction of test points...\n');
    
    %-- NEW
    model.mini = mini;
    %---
    [x, varx, modelUpdated] = vargplvmOptimisePoint(model, vardistx, Yts, display, iters);
    barmu = x;
    lambda = varx;
    [x, varx, modelUpdated] = vargplvmDynamicsUpdateModelTestVar(modelUpdated, barmu, lambda, Yts);
    
    % Find the absolute error
    Testmeans = x;
    Testcovars = varx;
    Varmu = vargplvmPosteriorMeanVar(modelUpdated, x, varx);
    % Mean error per pixel
    errorFull = sum(sum( abs(Varmu(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
    fprintf(1,'# GPLVM Error (in the missing dims) with missing inputs:%d\n', errorFull);
    
    errOnlyTimes2 = sum(sum( abs(Varmu2(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
    fprintf(1,'# GPLVM Error (in the missing dims) only times:%d\n', errOnlyTimes2);
    
    errorFullPr = sum(sum( abs(Varmu(:,indexPresent) - YtsOriginal(:,indexPresent)) ))/prod(size(YtsOriginal(:,indexPresent)));
    fprintf(1,'# GPLVM Error (in the present dims) with missing inputs:%d\n', errorFullPr);
    
    errorFullPr2 = sum(sum( abs(Varmu2(:,indexPresent) - YtsOriginal(:,indexPresent)) ))/prod(size(YtsOriginal(:,indexPresent)));
    fprintf(1,'# GPLVM Error (in the present dims) only times:%d\n', errorFullPr2);
    %     % Visualization of the reconstruction
    %     for i=1:Nstar
    %         subplot(1,2,1);
    %         fr=reshape(YtsOriginal(i,:),height,width);
    %         imagesc(fr);
    %         colormap('gray');
    %         subplot(1,2,2);
    %         fr=reshape(Varmu(i,:),height,width);
    %         imagesc(fr);
    %         colormap('gray');
    %         pause(0.5);
    %     end
    
    prunedModelUpdated = vargplvmPruneModel(modelUpdated,1);
    save([fileToSave(1:end-4) 'Pred.mat'], 'Testmeans', 'Testcovars', 'prunedModelUpdated');
    fprintf(1,'# Saved %s\n',[fileToSave(1:end-4) 'Pred.mat']);
    
    
    %-------- NN  ----------
    fprintf(1,'# NN(2) prediction...\n');
    mini =[];
    sortedInd = zeros(size(Yts,1), size(Ytr,1));
    for i=1:size(Yts,1)
        dst = dist2(Yts(i,indexPresent), Ytr(:,indexPresent));
        % mind: smaller distance %mini(i): smaller index
        %[mind, mini(i)] = min(dst);
        [void ,sortedInd(i,:)] = sort(dst,2);
    end
    clear dst %%%
    
    
    NNmuPart = zeros(Nstar, size(Ytr,2));
    % Set to 1 to find the two NN's in time space. set to 0 to find in
    % data space. timeNN=1 should be set only for datasets created with
    % the "everyTwo" option for the split.
    timeNN=0;
    k=2; % the k parameter in k-NN
    for i=1:Nstar
        if timeNN
            if i < Nstar
                NNmuPart(i, indexMissing) = 0.5*(Ytr(i,indexMissing) + Ytr(i+1,indexMissing));
            else
                NNmuPart(i, indexMissing) = Ytr(end,indexMissing);
            end
        else
            NNmuPart(i,indexMissing) = Ytr(sortedInd(i,1),indexMissing);
            for n=2:k
                NNmuPart(i,indexMissing) = NNmuPart(i,indexMissing)+Ytr(sortedInd(i,n),indexMissing);
            end
            NNmuPart(i,indexMissing) = NNmuPart(i,indexMissing)*(1/k);
        end
        NNmuPart(i, indexPresent) = Yts(i, indexPresent);
    end
    
    % Mean absolute error per pixel
    %errorNNPart = mean(abs(NNmuPart(:) - YtsOriginal(:)));
    errorNNPart = sum(sum( abs(NNmuPart(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
    fprintf(1,'# NN(2) Error (in the missing dims) with missing inputs:%d\n', errorNNPart);
    
    
    %%%%%%  Try different k for k-NN to see which is best
    fprintf(1,'# k-NN prediction for k=[1 3 4 5]...\n');
    if ~timeNN
        bestK=2; bestErr = errorNNPart; NNmuPartBest = NNmuPart;
        for k=[1 3 4 5]; % the k parameter in k-NN, try different values
            NNmuPartK = zeros(Nstar, size(Ytr,2));
            for i=1:Nstar
                NNmuPartK(i,indexMissing) = Ytr(sortedInd(i,1),indexMissing); % first NN
                for n=2:k % more NN's, if k>1
                    NNmuPartK(i,indexMissing) = NNmuPartK(i,indexMissing)+Ytr(sortedInd(i,n),indexMissing);
                end
                NNmuPartK(i,indexMissing) = NNmuPartK(i,indexMissing)*(1/k); % normalise with the number of NN's
                
                NNmuPartK(i, indexPresent) = Yts(i, indexPresent);
            end
            errorNNPartK = sum(sum( abs(NNmuPartK(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
            if errorNNPartK < bestErr
                bestErr = errorNNPartK;
                bestK=k;
                NNmuPartBest = NNmuPartK;
            end
        end
        clear NNmuPartK
    end
    % Mean absolute error per pixel
    fprintf(1,'# NNbest(%d) Error (in the missing dims) with missing inputs:%d\n',bestK, bestErr);
    %%%%%%%%%%%%%
    
    
    %     for i=1:Nstar
    %         subplot(1,2,1);
    %         fr=reshape(YtsOriginal(i,:),height,width);
    %         imagesc(fr);
    %         colormap('gray');
    %         subplot(1,2,2);
    %         fr=reshape(NNmuPart(i,:),height,width);
    %         imagesc(fr);
    %         colormap('gray');
    %         pause(0.5);
    %     end
end


%%
%--------------- SAMPLING -----%
if doSampling
    fprintf(1,'# Sampling...\n');
    if noiselessSamples
        [ySamp, xSamp] = vargpTimeDynamicsSample(model, noiselessSamples);
    else
        [ySamp, xSamp] = vargpTimeDynamicsSample(model);
    end
    bar(model.kern.comp{1}.inputScales)
    figure
    fprintf(1,'# Samples for X (press any key after each frame)...\n');
    for q=1:model.q
        % figure
        plot(model.dynamics.t_star, xSamp(:,q),'b-');
        hold on;
        plot(model.dynamics.t, model.X(:,q),'r-')
        hold off;
        legend('samples','original')
        title(['q=' num2str(q)])
        pause;
    end
    
    if ~exist('futurePred')
        errorSampling = mean(abs(ySamp(:) - YtsOriginal(:)));
        fprintf(1,'# Sampling Error (in the whole frame):%d\n', errorSampling);
        
    end
    
    %     fprintf('# Press any key to see the visualisation...\n');
    %     figure
    %     pause
    %     for i=1:Nstar
    %         subplot(1,2,1);
    %         fr=reshape(Varmu2(i,:),height,width);
    %         imagesc(fr);
    %         colormap('gray');
    %         subplot(1,2,2);
    %         fr=reshape(ySamp(i,:),height,width);
    %         imagesc(fr);
    %         colormap('gray');
    %         pause;
    %     end
end



%% Plots
plotLatentSpace=0;
if plotLatentSpace
    % find the largest dimensions
    [max, imax] = sort(model.kern.comp{1}.inputScales,'descend');
    
    N = size(model.X,1);
    Nstar = size(Yts,1);
    figure;
    fprintf(1,'The two larger latent dims after re-training with partial test data. Red are the test points\n');
    plot(model.X(1:model.N,imax(1)), model.X(1:model.N,imax(2)),'--rs', 'Color','b');
    hold on
    if exist('Testmeans2')
        plot(Testmeans2(:,imax(1)),Testmeans2(:,imax(2)),'--rs', 'Color','r');
    end
    if exist('xSamp')
        plot(xSamp(:,imax(1)),xSamp(:,imax(2)),'--rs', 'Color','g');
    end
    title('Visualization of the latent space after re-training with partial test data');
    hold off
end


%%
if ~exist('resultsDynamic')  resultsDynamic = 0; end
if resultsDynamic
    % The following causes OUTOFMEMORY exception except from when we prune the
    % video dimensions: (Note: also, lvmVisualise was changed a bit so that no
    % dynamic slides is presented, because otherwise a strange error occurs).
    model.y = Ytr;
    [modelP, newHeight, newWidth] = vargplvmReduceVidModel(model, height, width, 4,4);
    lvmVisualise(modelP, [], 'imageVisualise', 'imageModify', [newHeight newWidth],1,0,1);
    clear modelP
end
