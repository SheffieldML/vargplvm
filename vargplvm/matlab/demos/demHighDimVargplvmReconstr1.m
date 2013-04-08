% This script is no longer used. See demHighDimVargplvmTrained.m instead
% parameters.

% TODO: make the mask not to be the half of the image but randomly selected
% pixels (50% missing).


%----- Constants
randn('seed', 1e5);
rand('seed', 1e5);

% load data
dataSetName = 'missa'; % 360x288
%load miss-americaHD
Y = lvmLoadData(dataSetName);

width=360;
height=288;


if ~exist('trainedModel')      trainedModel = 'demMissaVargplvm2';      end
if ~exist('reconstrIters')     reconstrIters = 100;      end
displayMovies = 0;

%------- LOAD THE MODEL


N = size(Y,1);
% assumes the number of data is even
Nstar = N/2; Nts = N/2; Ntr = N/2;
indTr = 1:2:N; indTs = 2:2:N;
Ytr = Y(indTr,:); Yts = Y(indTs,:); YtsOriginal = Yts;
t = linspace(0, 2*pi, size(Y, 1)+1)'; t = t(1:end-1, 1);
timeStampsTraining = t(indTr,1); timeStampsTest = t(indTs,1);

load(trainedModel)
if prunedModel.fixInducing
    prunedModel.inducingIndices=1:size(Ytr,1);
    prunedModelTr.inducingIndices=1:size(Ytr,1);
    prunedModelInit.inducingIndices=1:size(Ytr,1);
end

model = vargplvmRestorePrunedModel(prunedModel, Ytr);
modelInit = vargplvmRestorePrunedModel(prunedModelInit, Ytr);
modelTr = vargplvmRestorePrunedModel(prunedModelTr, Ytr);

model = modelTr;
model.y = Ytr;
model.m= gpComputeM(model); %%%
model.y=[];


%%

%--------------------- Prediction using the only information in the test time points

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

fprintf(1,'# Error GPLVM only times: %d\n', errorOnlyTimes);
fprintf(1,'# Error NN (no partial inf.): %d\n', errorNN);


%%

%-------------  Reconstruction With Partial Observations ---------------
%--- Remove some dimensions from the test points

display = 1;
indexP = [];
Init = [];


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

%%

% See the result
if displayMovies
    for i=1:size(Ytemp,1)
        fr=reshape(Ytemp(i,:),288,360);
        imagesc(fr);
        colormap('gray');
        pause(0.1)
    end
    clear Ytemp
end
%%

%------- Predict with vargplvm
mini =[];
for i=1:size(Yts,1)
    % initialize the latent points using the nearest neighbour
    % from he training data
    dst = dist2(Yts(i,indexPresent), Ytr(:,indexPresent));
    [mind, mini(i)] = min(dst);
end

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
fprintf(1,'# Vargplvm Error (in the missing dims) with missing inputs:%d\n', errorFull);

%%

% Visualization of the reconstruction
if displayMovies
    for i=1:Nstar
        subplot(1,2,1);
        fr=reshape(YtsOriginal(i,:),288,360);
        imagesc(fr);
        colormap('gray');
        subplot(1,2,2);
        fr=reshape(Varmu(i,:),288,360);
        imagesc(fr);
        colormap('gray');
        pause(0.5);
    end
end

%%
%---------- Prediction with NN

% The following code finds the NN in the data space. However, if we print
% dst we see that tne NN in the data space is almost the same as the NN in
% the time space, i.e. the next or previous frame.
% mini =[];
% for i=1:size(Yts,1)
%     dst = dist2(Yts(i,indexPresent), Ytr(:,indexPresent));
%     [mind, mini(i)] = min(dst);
% end
%

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

%%

% Visualization of the reconstruction
if displayMovies
    for i=1:Nstar
        subplot(1,2,1);
        fr=reshape(YtsOriginal(i,:),288,360);
        imagesc(fr);
        colormap('gray');
        subplot(1,2,2);
        fr=reshape(NNmuPart(i,:),288,360);
        imagesc(fr);
        colormap('gray');
        pause(0.5);
    end
end
%%
