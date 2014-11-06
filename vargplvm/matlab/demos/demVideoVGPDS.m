randn('seed', 1e6);
rand('seed', 1e6);
if ~exist('mappingKern', 'var'), mappingKern = 'rbfardjit'; end
if ~exist('indPoints', 'var'), indPoints = 20; end
if ~exist('latentDim', 'var'), latentDim = 10; end
% Define a temporal model
if ~exist('dynamicsConstrainType', 'var'), dynamicsConstrainType = {'time'}; end
if ~exist('initVardistIters', 'var'), initVardistIters = 500; end
if ~exist('itNo', 'var'), itNo = [50 50]; end
if ~exist('dynamicKern', 'var'), dynamicKern = {'matern32','white','bias'}; end


%%
% load data
[Y, lbls] = lvmLoadData('brendan');
Y = Y(801:880,:);
t = linspace(1,2*pi, size(Y,1))';

h = 20;
w = 28;
scrsz = get(0,'ScreenSize');
figure('Position',[0.2*scrsz(3) 0.2*scrsz(4) 20*5 28*5])
for i=1:size(Y,1)
    imagesc(reshape(Y(i,:), h, w)'); colormap('gray');
    pause(0.05);
end

%% Split between training and test blocks
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
Nstar = length(indTs);

Ytr = Y(indTr,:); t_tr = t(indTr, :);
Yts = Y(indTs,:); t_ts = t(indTs, :);

% Additionally define some pixels to be observed for the test set
cutPoint=round(h/2)+1;
mask = zeros(h,1);
mask(cutPoint) = 1;
mask = [ones(1,cutPoint) zeros(1,h-cutPoint)];
mask=repmat(mask, 1,w);
indexMissing = find(mask);

YtsOriginal = Yts;
Yts(:, indexMissing) = NaN;
% Play the test set fully observed (left) and partially observed (right)
scrsz = get(0,'ScreenSize');
figure('Position',[0.2*scrsz(3) 0.2*scrsz(4) 20*15 28*5])
for i=1:size(Yts,1)
    subplot(1,2,1)
    imagesc(reshape(YtsOriginal(i,:), h, w)'); colormap('gray');
    subplot(1,2,2);
    imagesc(reshape(Yts(i,:), h, w)'); colormap('gray');
    pause(0.05);
end
%% Train models

% Options for the model
vargplvm_init; % Returns a configuration structure 'globalOpt'
options = vargplvmOptions('dtcvar');
options.kern = mappingKern;
options.numActive = indPoints;
options.optimiser = 'scg2';
options.latentDim = latentDim;
options.initSNR = 100;
if ~isempty('dynamicsConstrainType')
    % Temporal model (VGPDS)
    optionsDyn.type = 'vargpTime';
    optionsDyn.t=t;
    if ~isstruct(dynamicKern)
        optionsDyn.kern = kernCreate(t, dynamicKern);
    else
        optionsDyn.kern = dynamicKern;
    end
    % Create and optimise the model
    [~, ~, ~, ~, modelInitVardist] = vargplvmEmbed2(Ytr, latentDim, options, initVardistIters, 0, true, optionsDyn);
     model = vargplvmOptimiseModel(modelInitVardist, true, true, {0,itNo}, true);
else
    % Non-dynamical model (Bayesian GP-LVM)
    % Create and optimise the model
    [~, ~, ~, ~, modelInitVardist] = vargplvmEmbed2(Ytr, latentDim, options, initVardistIters, 0, true);
    model = vargplvmOptimiseModel(modelInitVardist, true, true, {0,itNo}, true);
end

%%

fprintf(1, '# Predicting only with the test time points...\n')
[Testmeans2 Testcovars2] = vargplvmPredictPoint(model.dynamics, t_ts);
Varmu2 = vargplvmPosteriorMeanVar(model, Testmeans2, Testcovars2);
% Mean absolute error per pixel
errorOnlyTimes = mean(abs(Varmu2(:) - YtsOriginal(:)));

%%
% Play the test set fully observed (left) and predicted (right)
scrsz = get(0,'ScreenSize');
figure('Position',[0.2*scrsz(3) 0.2*scrsz(4) 20*15 28*5])
t_ind = 1:size(Y,1);
t_ind = t_ind(indTs);
for i=1:size(Yts,1)
    subplot(1,2,1)
    imagesc(reshape(YtsOriginal(i,:), h, w)'); colormap('gray'); title(num2str(t_ind(i)))
    subplot(1,2,2);
    imagesc(reshape(Varmu2(i,:), h, w)'); colormap('gray');
    pause;
end
