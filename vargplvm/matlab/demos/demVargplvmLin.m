


% DEMOIL100VARGPLVM1 Run variational GPLVM on 100 points from the oil data.

% SHEFFIELDML
clear
close all

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

experimentNo = 1;
printDiagram = 1;

initVardistIters = 200;
iters = 200;
Q=12;
K=25;

dataSetName = 'oil100';
[Y,lbls] = vargplvmLoadData(dataSetName);
%Y = Y(1:20,:);
% dataSetName = 'bc/wine';
% addpath(genpath('../../bc-vargplvm/matlab'));
% globalOpt.dataSetName = dataSetName;
% globalOpt.dataPerClass = 30;
% globalOpt.dataToKeep = -1;
% [globalOpt, Y,lbls] = bc_LoadData(globalOpt);
% lbls = Y.lbls;
% Y = Y.Y;
% dataSetName = 'wine';

% Set up model
options = vargplvmOptions('dtcvar');
options.kern = {'linard2','white'}; %{'linard2', 'bias','white'}; % 'rbfardjit'; %{'rbfard2','white','bias'};%'rbfardjit';
options.numActive = K;
Xpca = ppcaEmbed(Y,Q);
options.initX = Xpca; %rand(size(Xpca));
%options.tieParam = 'tied';

options.optimiser = 'scg';%'scg2';
latentDim = Q;
d = size(Y, 2);

% demo using the variational inference method for the gplvm model
model = vargplvmCreate(latentDim, d, Y, options);
%
%model = vargplvmParamInit(model, model.m, model.X);
model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));
model.beta = 100 / var(Y(:));
%model.vardist.covars = 0.2*ones(size(model.vardist.covars));

%--- TEMP
model.kern.comp{1}.inputScales = model.kern.comp{1}.inputScales./model.kern.comp{1}.inputScales;

% Optimise the model.
display = 1;
%model.fixParamIndices = [2301:2310];
modelInit = model;

%%

model.learnBeta = 0; model.learnSigmaf=0;
model = vargplvmOptimise(model, display, initVardistIters);
model.learnBeta = 1; model.learnSigmaf=1;

model = vargplvmOptimise(model, display, iters);
fprintf('1/model.beta = %.5f, var(Y(:))=%.5f\n', 1/model.beta, var(Y(:)));

%%
capName = dataSetName;;
capName(1) = upper(capName(1));
modelType = model.type;
modelType(1) = upper(modelType(1));
save(['dem' capName modelType num2str(experimentNo) '.mat'], 'model');

% order wrt to the inputScales
mm = vargplvmReduceModel(model,2);
modelInit.kern.comp{1}.inputScales = model.kern.comp{1}.inputScales;
mmInit = vargplvmReduceModel(modelInit,2);
% plot the two largest twe latent dimensions
if exist('printDiagram') & printDiagram
    % Plot Vargplvm:
    vargplvmPrintPlot(mm, lbls, capName, experimentNo,false);
    title('vargplvm')
    % PLot PCA:
    mm.XOrig = mm.X;
    mm.X = ppcaEmbed(Y,2);
    vargplvmPrintPlot(mm, lbls, capName, experimentNo,false);
    title('pca')
    % VargplvmInit
    vargplvmPrintPlot(mmInit, lbls, capName, experimentNo,false);
    title('INIT vargplvm')
    % Plot scales
    figure
    subplot(1,2,1)
    bar(sort(model.kern.comp{1}.inputScales, 'descend'))
    subplot(1,2,2)
    bar(pca(Y,10))
    
    mm.X = mm.XOrig;
end
% For 3D plots:
%labels = transformLabels(lbls); dims = [1 2 3];
%plot3k({model.X(:,dims(1)) model.X(:,dims(2)) model.X(:,dims(3))}, 'ColorData', labels, 'Marker', {'x',6});

