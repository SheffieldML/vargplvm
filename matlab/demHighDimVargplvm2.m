% DEMHIGHDIMVARGPLVM2
% This is an old script. See demHighDimVargplvm3.m
% VARGPLVM
clear timeStamps; % in case it's left from a previous experiment

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

% Define constants (in a manner that allows other scripts to parametrize
% this one).
if ~exist('experimentNo')   experimentNo = 404;      end
if ~exist('itNo')           itNo =2000;              end     % Default: 2000
if ~exist('indPoints')      indPoints = 49;          end     % Default: 49
if ~exist('latentDim')      latentDim = 40;          end
% Set to 1 to use dynamics or to 0 to use the standard var-GPLVM
if ~exist('dynUsed')        dynUsed = 1;             end
% Set to 1 to keep only the dimensions modelling the head (missa dataset)
if ~exist('cropVideo')      cropVideo = 0;           end
% Set to 1 to remove the frames with the translation (missa dataset)
if ~exist('removeTransl')   removeTransl = 0;        end
if ~exist('fixedBetaIters') fixedBetaIters = 0;      end     % DEFAULT: 23
% Set to 1 to tie the inducing points with the latent vars. X
if ~exist('fixInd') fixInd = 0;                      end
if ~exist('dynamicKern') dynamicKern = {'rbf', 'white', 'bias'}; end
if ~exist('reconstrIters') reconstrIters = 1000;                 end
if ~exist('vardistCovarsMult') vardistCovarsMult=1.3;            end % 0.1 gives around 0.5 init.covars. 1.3 biases towards 0.

        
% load data
dataSetName = 'ocean'; % 720x1280

load '../../VideoManipulation/datasets/mat/DATA_Ocean';
% For this dataset there is a translation in space between frames 66-103


% dataSetName = 'susie';
% load susie_352_240

fprintf(1,'\n#----------------------------------------------------\n');
fprintf(1,'# Dataset: %s\n',dataSetName);
fprintf(1,'# ExperimentNo: %d\n', experimentNo);
fprintf(1,'# Inducing points: %d\n',indPoints);
fprintf(1,'# Latent dimensions: %d\n',latentDim);
fprintf(1,'# Iterations (with/without fixed Beta): %d / %s\n',fixedBetaIters, num2str(itNo));
fprintf(1,'# Tie Inducing points: %d\n',fixInd);
fprintf(1,'# Dynamics used: %d\n', dynUsed);
if dynUsed
    fprintf(1,'# Dynamics kern: ');
    disp(dynamicKern);
end

fprintf(1,'# CropVideo / removeTranslation: %d / %d \n', cropVideo, removeTransl);
fprintf(1,'# VardistCovarsMult: %d \n', vardistCovarsMult);
%fprintf(1,'#----------------------------------------------------\n');

switch dataSetName
    case 'missa'
        width=360;
        height=288;
    case 'susie'
        width=352;
        height=240;
    case 'ocean'
        width=1280;
        height=720;
end

% Remove the translation part (only for the "missa" dataset)
if removeTransl
    Y = [Y(1:64,:) ; Y(103:end,:)];
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
    end
end


% %%%%TEMP (for gradchek)
% Y = Y(1:30,:);
% %%%%%

dims = size(Y,2);

trainModel = 1;% Set it to 1 to retrain the model. Set it to 0 to load an already trained one.

% Take a downsample version of the video for training and test on the
% remaining frames
N = size(Y,1);
% assumes the number of data is even
Nstar = N/2; Nts = N/2; Ntr = N/2;
indTr = 1:2:N; indTs = 2:2:N;
Ytr = Y(indTr,:); Yts = Y(indTs,:); YtsOriginal = Yts;
t = linspace(0, 2*pi, size(Y, 1)+1)'; t = t(1:end-1, 1);
timeStampsTraining = t(indTr,1); timeStampsTest = t(indTs,1);

fprintf(1,'# Dataset size used (train/test) : %d / %d \n', size(Ytr,1), size(Yts,1));
%fprintf(1,'# CropVideo / removeTranslation: %d / %d \n', cropVideo, removeTransl);
fprintf(1,'#----------------------------------------------------\n');


%%
% % Play movie
% for i=1:size(Y,1)
%     fr=reshape(Y(i,:),height,width);
%     imagesc(fr); colormap('gray');
%     pause(0.08);
% end

%%


% Set up model
options = vargplvmOptions('dtcvar');
options.kern = {'rbfard2', 'bias', 'white'};
options.numActive = indPoints;
options.optimiser = 'scg';
d = size(Y, 2);


if trainModel
    % demo using the variational inference method for the gplvm model
    fprintf(1,'# Creating the model...\n');
    % scale = std(Ytr);
    % scale(find(scale==0)) = 1;
    %options.scaleVal = mean(std(Ytr));
    options.scaleVal = sqrt(var(Ytr(:)));
    if fixInd
        options.fixInducing=1;
        options.fixIndices=1:size(Ytr,1);
    end
    model = vargplvmCreate(latentDim, d, Ytr, options);
    
    
    
    % Temporary: in this demo there should always exist the mOrig field
    if ~isfield(model, 'mOrig')
        model.mOrig = model.m;
    end
    
    model = vargplvmParamInit(model, model.mOrig, model.X);
    model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));
    %model.kern.comp{1}.variance = max(var(Y)); %%%
    
    
    %-------- Add dynamics to the model -----
    if dynUsed
        fprintf(1,'# Adding dynamics to the model...\n');
        optionsDyn.type = 'vargpTime';
        optionsDyn.t=timeStampsTraining;
        optionsDyn.inverseWidth=100; % Default: 100
                
        
        kern = kernCreate(t, dynamicKern); % Default: {'rbf','white','bias'}
        kern.comp{2}.variance = 1e-1; % Usual values: 1e-1, 1e-3
        % The following is related to the expected number of
        % zero-crossings.(larger inv.width numerator, rougher func)
        if ~strcmp(kern.comp{1}.type,'ou')
            kern.comp{1}.inverseWidth = optionsDyn.inverseWidth./(((max(t)-min(t))).^2);
            kern.comp{1}.variance = 1;
        end
        optionsDyn.kern = kern;
        optionsDyn.vardistCovars = vardistCovarsMult; % 0.23 gives true vardist.covars around 0.5 (DEFAULT: 0.23) for the ocean dataset
             
        % Fill in with default values whatever is not already set
        optionsDyn = vargplvmOptionsDyn(optionsDyn);
        model = vargplvmAddDynamics(model, 'vargpTime', optionsDyn, optionsDyn.t, 0, 0,optionsDyn.seq);
        
        fprintf(1,'# Further calibration of the initial values...\n');
        model = vargplvmInitDynamics(model,optionsDyn);
    end
    
    model.beta=1/(0.01*var(model.mOrig(:)));
    modelInit = model;
    
    %%---
    capName = dataSetName;
    capName(1) = upper(capName(1));
    modelType = model.type;
    modelType(1) = upper(modelType(1));
    fileToSave = ['dem' capName modelType num2str(experimentNo) '.mat'];
    %%---
    
    display = 1;
    %     %%%% Optimisation
    %     % do not learn beta for few iterations for intitialization
    if fixedBetaIters ~=0
        model.learnBeta = 0;
        fprintf(1,'# Intitiliazing the model (fixed beta)...\n');
        model = vargplvmOptimise(model, display, fixedBetaIters); % Default: 20
        fprintf(1,'1/b = %.4d\n',1/model.beta);
        
        %fprintf(1,'# Saving %s\n',fileToSave);
        %prunedModel = vargplvmPruneModel(model);
        %prunedModelInit = vargplvmPruneModel(modelInit);
        %save(fileToSave, 'prunedModel', 'prunedModelInit');
        
        modelBetaFixed = model;
        model.fixedBetaIters=fixedBetaIters;
    end
    
    model.learnBeta = 1;
    model.iters = 0;
    
    % Optimise the model.
    for i=1:length(itNo)
        iters = itNo(i); % default: 2000
        fprintf(1,'\n# Optimising the model for %d iterations (session %d)...\n',iters,i);
        model = vargplvmOptimise(model, display, iters);
        model.iters = model.iters + iters;
        fprintf(1,'1/b = %.4d\n',1/model.beta);
        modelTr = model;
        
        fprintf(1,'# 1/b=%.4f\n var(m)=%.4f\n 1/init_b=%.4f\n',1/model.beta, var(model.mOrig(:)), modelInit.beta);
        
        % Save model
        fprintf(1,'# Saving %s\n',fileToSave);
        prunedModel = vargplvmPruneModel(model);
        prunedModelInit = vargplvmPruneModel(modelInit);
        prunedModelTr = vargplvmPruneModel(modelTr);
        save(fileToSave, 'prunedModel', 'prunedModelInit', 'prunedModelTr');
    end
else
    % Load pre-trained model
    capName = dataSetName;
    capName(1) = upper(capName(1));
    modelType = 'vargplvm';
    modelType(1) = upper(modelType(1));
    fileToLoad = ['dem' capName modelType num2str(experimentNo) '.mat'];
    load(fileToLoad);
    fileToSave = fileToLoad;
end


% Don't make predictions for the non-dynamic case
if ~dynUsed
    return
end



%%%-----------------------   RECONSTRUCTION --------------------%%%

model = modelTr;
model.y = Ytr;
model.m= gpComputeM(model); %%%
model.y=[];

% Prediction using the only information in the test time points
[Testmeans2 Testcovars2] = vargplvmPredictPoint(model.dynamics, timeStampsTest);
Varmu2 = vargplvmPosteriorMeanVar(model, Testmeans2, Testcovars2);
% Mean absolute error per pixel
errorOnlyTimes = mean(abs(Varmu2(:) - YtsOriginal(:)));

% 2-nearest (in time) prediction
NNmu = zeros(Nstar, model.d);
for i=1:Nstar
    if i < Nstar
        NNmu(i, :) = 0.5*(Ytr(i,:) + Ytr(i+1,:));
    else
        NNmu(i, :) = Ytr(end,:);
    end
end
% Mean absolute error per pixel
errorNN = mean(abs(NNmu(:) - YtsOriginal(:)));
%errorNN = sum(sum( abs(NNmu - YtsOriginal)) )/prod(size(YtsOriginal)); %equivalent

fprintf(1,'# Error GPLVM: %d\n', errorOnlyTimes);
fprintf(1,'# Error NN: %d\n', errorNN);

return

% Visualization of the reconstruction
for i=1:Nstar
    subplot(1,2,1);
    fr=reshape(YtsOriginal(i,:),height,width);
    imagesc(fr);
    colormap('gray');
    subplot(1,2,2);
    fr=reshape(Varmu2(i,:),height,width);
    imagesc(fr);
    colormap('gray');
    pause(0.2);
end





%%%------------------ Extra --------------------%%%%%%%

% The full movie!!
for i=1:size(Varmu2,1)
    fr=reshape(Ytr(i,:),height,width);
    imagesc(fr);
    colormap('gray');
    pause(0.09)
    fr=reshape(Varmu2(i,:),height,width);
    imagesc(fr);
    colormap('gray');
    pause(0.001)
end
%%%%%





% if you like you can also do prediction after training with partially
% observed test images (also good for de-bugging because this should be better
% than the onlyTimes prediction)

predWithMs = 0;

if predWithMs
    %
    display = 1;
    indexP = [];
    Init = [];
    
    fprintf(1, '# Partial reconstruction of test points...\n');
    w=360; % width
    h=288; % height
    cutPoint=h/2+13;
    mask = [ones(1,cutPoint) zeros(1,h-cutPoint)];
    mask=repmat(mask, 1,w);
    indexMissing = find(mask);
    indexPresent = setdiff(1:model.d, indexMissing);
    
    %
    Yts = YtsOriginal;
    Ytemp=Yts;
    Ytemp(:,indexMissing)=NaN;
    for i=1:size(Ytemp,1)
        fr=reshape(Ytemp(i,:),288,360);
        imagesc(fr);
        colormap('gray');
        pause(0.1)
    end
    clear Ytemp
    
    mini =[];
    for i=1:size(Yts,1)
        % initialize the latent points using the nearest neighbour
        % from he training data
        dst = dist2(Yts(i,indexPresent), Ytr(:,indexPresent));
        [mind, mini(i)] = min(dst);
    end
    
    % new barmu
    %jit = 1e-6;
    %%muInit = [model.vardist.means; model.vardist.means(mini,:)];
    %Kt = kernCompute(model.dynamics.kern, timeStampsTest);
    %Lkt = chol(Kt + jit*mean(diag(Kt))*eye(size(Kt,1)))';
    %%barmuInit = Lkt'\(Lkt\muInit);
    %%model.dynamics.vardist.means = barmuInit(1:model.N,:);
    %vardistx = vardistCreate(model.vardist.means(mini,:), model.q, 'gaussian');
    %vardistx.means = Lkt'\(Lkt\model.vardist.means(mini,:)) + 0.5*randn(size(vardistx.means));
    %vardistx.covars = 0.2*ones(size(vardistx.covars));
    vardistx = vardistCreate(model.dynamics.vardist.means(mini,:), model.q, 'gaussian');
    vardistx.covars = model.dynamics.vardist.covars(mini,:);%0.2*ones(size(vardistx.covars));
    model.vardistx = vardistx;
    
    Yts(:,indexMissing) = NaN;
    model.dynamics.t_star = timeStampsTest;
    iters=reconstrIters;
    [x, varx, modelUpdated] = vargplvmOptimisePoint(model, vardistx, Yts, display, iters);
    barmu = x;
    lambda = varx;
    [x, varx, modelUpdated] = vargplvmDynamicsUpdateModelTestVar(modelUpdated, barmu, lambda, Yts);
    
    % Find the absoluate error
    Testmeans = x;
    Testcovars = varx;
    Varmu = vargplvmPosteriorMeanVar(modelUpdated, x, varx);
    % Mean error per pixel
    errorFull = sum(sum( abs(Varmu(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
    fprintf(1,'# GPLVM Error (in the missing dims) with missing inputs:%d\n', errorFull);
    
    % Visualization of the reconstruction
    for i=1:Nstar
        subplot(1,2,1);
        fr=reshape(YtsOriginal(i,:),height,width);
        imagesc(fr);
        colormap('gray');
        subplot(1,2,2);
        fr=reshape(Varmu(i,:),height,width);
        imagesc(fr);
        colormap('gray');
        pause(0.5);
    end
    
    NNmuPart = zeros(Nstar, model.d);
    for i=1:Nstar
        if i < Nstar
            NNmuPart(i, indexMissing) = 0.5*(Ytr(i,indexMissing) + Ytr(i+1,indexMissing));
        else
            NNmuPart(i, indexMissing) = Ytr(end,indexMissing);
        end
        NNmuPart(i, indexPresent) = Yts(i, indexPresent);
    end
    % Mean absolute error per pixel
    %errorNNPart = mean(abs(NNmuPart(:) - YtsOriginal(:)));
    errorNNPart = sum(sum( abs(NNmuPart(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
    fprintf(1,'# NN Error (in the missing dims) with missing inputs:%d\n', errorNNPart);
    
end
