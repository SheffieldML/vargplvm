
% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);


%--
dataSetName = 'weizmann_horses';
if ~exist('latentDim'), latentDim = 15; end
if ~exist('indPoints'), indPoints = 80; end
if ~exist('initVardistIters'), initVardistIters = 500; end
if ~exist('itNo'), itNo = 1000; end
if ~exist('scale2var1'), scale2var1 = true; end
%--

vargplvm_init;
if ~isempty(globalOpt.diaryFile)
    diary(globalOpt.diaryFile)
end

[Ytr, lbls, Yts] = vargplvmLoadData(globalOpt.dataSetName);
height = lbls(1); width = lbls(2);
YtsOriginal = Yts;

% Binarize values
zeroInd = find(Ytr < 255/2);
oneInd = find(Ytr >= 255/2);
Ytr(zeroInd) =  0;
Ytr(oneInd) = 1;
clear('zeroInd', 'oneInd');

% playMov(height, width, 0.08,Ytr)
% playMov(height, width, 0.12,Yts)

%%

% Set up model
options = vargplvmOptions('dtcvar');
options.kern = globalOpt.mappingKern;
options.numActive = globalOpt.indPoints;
options.optimiser = 'scg2';

options.enableDgtN = globalOpt.DgtN;



% scale = std(Ytr);
% scale(find(scale==0)) = 1;
%options.scaleVal = mean(std(Ytr));
options.scaleVal = sqrt(var(Ytr(:)));
if globalOpt.fixInd
    options.fixInducing=1;
    options.fixIndices=1:size(Ytr,1);
end
options.scale2var1 = globalOpt.scale2var1;
fprintf(1,'# Creating the model...\n');
model = vargplvmCreate(globalOpt.latentDim, size(Ytr,2), Ytr, options);
model = vargplvmModelInit(model, globalOpt);
modelInit = model; %%%


%%
model = vargplvmOptimiseModel(model, true, true, {globalOpt.initVardistIters,globalOpt.itNo});

%vargplvmWriteResult(model, model.type, globalOpt.dataSetName, globalOpt.experimentNo);
%if exist('printDiagram') & printDiagram
%    lvmPrintPlot(model, [], globalOpt.dataSetName, globalOpt.experimentNo);
%end

% Load the results and display dynamically.
%lvmVisualise(model,[],'imageVisualise','imageModify', [32 32], 1, 0, 1);

%%

%%%-----------------------   RECONSTRUCTION --------------------%%%

%%
if ~globalOpt.doPredictions
    if ~isempty(globalOpt.diaryFile)
        diary off
    end
    return
end

%clear modelTr %%
%model.y = Ytr;
%model.m= gpComputeM(model); %%%
%model.y=[];
Nstar = size(YtsOriginal,1);

%
display = 1;
indexP = [];
Init = [];

fprintf(1, '# Partial reconstruction of test points...\n');
% w=360; % width
% h=288; % height
w=width; h=height;

if ~exist('cut')
    cut = 'vertically';
end

switch cut
    case 'horizontally'
        cutPoint=round(h/2)+13;
        if exist('cutP')    cutPoint = cutP; end
        mask = [ones(1,cutPoint) zeros(1,h-cutPoint)];
        mask=repmat(mask, 1,w);
        indexMissing = find(mask);
    case 'vertically'
        indexMissing=1:round(size(Yts,2)/2);
end
indexPresent = setdiff(1:model.d, indexMissing);

%
Yts = YtsOriginal;
% See the dataset after cutting some dims
%{
    Ytemp=Yts;
    Ytemp(:,indexMissing)=NaN;
    for i=1:size(Ytemp,1)
        fr=reshape(Ytemp(i,:),h,w);
        imagesc(fr);
        colormap('gray');
        pause
    end
    clear Ytemp
%}


Yts(:,indexMissing) = NaN;

indexP = [];
Init = [];
Testmeans = [];
Testcovars = [];
Varmu = [];
Varsigma = [];

fprintf(1, '# Partial reconstruction of test points...\n');

% patrial reconstruction of test points
for i=1:size(Yts,1)
    %indexPresent = indicesPresent(i,:);
    indexP(i,:) = indexPresent;
    
    % initialize the latent point using the nearest neighbour
    % from he training data
    dst = dist2(Yts(i,indexPresent), Ytr(:,indexPresent));
    [mind, mini] = min(dst);
    
    Init(i,:) = model.vardist.means(mini,:);
    % create the variational distribtion for the test latent point
    vardistx = vardistCreate(model.vardist.means(mini,:), model.q, 'gaussian');
    vardistx.covars = 0.2*ones(size(vardistx.covars));
    
    % optimize mean and vars of the latent point
    model.vardistx = vardistx;
    [x, varx] = vargplvmOptimisePoint(model, vardistx, Yts(i, :), display, globalOpt.reconstrIters);
    
    Testmeans(i,:) = x;
    Testcovars(i,:) = varx;
    
    % reconstruct the missing outputs
    [mu, sigma] = vargplvmPosteriorMeanVar(model, x, varx);
    Varmu(i,:) = mu;
    Varsigma(i,:) = sigma;
    
    %
end


VarmuOrig = Varmu;
Varmu(Varmu>=0.4)=1;
Varmu(Varmu<0.4)=0;

% Mean error per pixel
errorFull = sum(sum( abs(Varmu(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
fprintf(1,'# GPLVM Error (in the missing dims) with missing inputs:%d\n', errorFull);

errorFullPr = sum(sum( abs(Varmu(:,indexPresent) - YtsOriginal(:,indexPresent)) ))/prod(size(YtsOriginal(:,indexPresent)));
fprintf(1,'# GPLVM Error (in the present dims) with missing inputs:%d\n', errorFullPr);



%-------- NN  ----------
fprintf(1,'# NN prediction...\n');
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
k=2; % the k parameter in k-NN
for i=1:Nstar
    NNmuPart(i,indexMissing) = Ytr(sortedInd(i,1),indexMissing);
    for n=2:k
        NNmuPart(i,indexMissing) = NNmuPart(i,indexMissing)+Ytr(sortedInd(i,n),indexMissing);
    end
    NNmuPart(i,indexMissing) = NNmuPart(i,indexMissing)*(1/k);
    NNmuPart(i, indexPresent) = Yts(i, indexPresent);
end

% Mean absolute error per pixel
%errorNNPart = mean(abs(NNmuPart(:) - YtsOriginal(:)));
errorNNPart = sum(sum( abs(NNmuPart(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
fprintf(1,'# NN(2) Error (in the missing dims) with missing inputs:%d\n', errorNNPart);


%%%%%%  Try different k for k-NN to see which is best

bestK=2; bestErr = errorNNPart; NNmuPartBest = NNmuPart;
clear NNmuPart %%%
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
% Mean absolute error per pixel
fprintf(1,'# NNbest(%d) Error (in the missing dims) with missing inputs:%d\n',bestK, bestErr);
%%%%%%%%%%%%%

Varmu2 = VarmuOrig;
Varmu2(Varmu2>0.75)=1;
Varmu2(Varmu2<0.25)=0;
% Means
for i=1:size(Yts,1)
    subplot(1,5,1)
    imagesc(reshape(VarmuOrig(i,:), 32,32)), colormap('gray'), title('GP-LVM')
    subplot(1,5,2)
    imagesc(reshape(Varmu2(i,:), 32,32)), colormap('gray'), title('GP-LVM')
    subplot(1,5,3)
    imagesc(reshape(Ytr(sortedInd(i,1),:), 32,32)), colormap('gray'), title('NN')
    subplot(1,5,4)
    imagesc(reshape(Yts(i,:),32,32)), colormap('gray'), title('Test')
    subplot(1,5,5)
    imagesc(reshape(YtsOriginal(i,:),32,32)), colormap('gray'), title('Original')
    pause
end


% Samples
numSamples = 5;
for i=1:size(Yts,1)
    subplot(1,numSamples+1,1)
    imagesc(reshape(YtsOriginal(i,:), height, width)), colormap('gray'), title('original')
    for j=1:numSamples
        xSamp = Testmeans(i,:) + sqrt(Testcovars(i,:)).*randn(size(Testmeans(i,:)));
        %xSamp = gsamp(zeros(1, size(diag(Testcovars(i,:)), 1)), diag(Testcovars(i,:)), 1);
        %xSamp(3:end) = 0;
        [muSamp, varSamp] = vargplvmPosteriorMeanVar(model, xSamp);
        ySamp = randn(size(muSamp)).*sqrt(varSamp)+ muSamp;
        %ySamp(ySamp  < 0.1)=0; ySamp(ySamp > 0.8)=1; %%
        subplot(1,numSamples+1,j+1)
        imagesc(reshape(muSamp, height, width)), colormap('gray')
    end
    hold on
    pause
end


if ~isempty(globalOpt.diaryFile)
    diary off
end