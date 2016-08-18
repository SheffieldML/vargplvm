%%%%%%%%%%%%%%%%%% HOW TO RUN THIS DEMO %%%%%%%%%%%%%
% Invoke it with the following commands / suggested parameters:

%{
 experimentNo = 1;
 indEvery = 20; % Place one inducing point every 20 test points
 myRandSeed = 1001;
 dataSetName = 'mgdata';
 windowSize = 18; % Window size for autoregressive prediction
 windowsTr = 4;   % How many windondows of training data to use
 windowsExtrap = ceil(500/windowSize);  % How many windows to predict when extrapolating
 parallel = true; % If you have the parallel toolbox
 evaluation = true;
 demKstepAheadMg_master;
%}
%%%%%%%%%%%%%%%%%%%%%%%%



if ~exist('myRandSeed', 'var'), myRandSeed = floor(rand*1e5); end
randn('seed', myRandSeed); rand('seed', myRandSeed);


if ~exist('experimentNo', 'var'), experimentNo = 404; end
if ~exist('dataSetName', 'var'), dataSetName = 'mgdata'; end
if ~exist('windowsTr', 'var'), windowsTr = 18; end
if ~exist('windowsInterp', 'var'), windowsInterp = 0; end
if ~exist('windowsExtrap', 'var'), windowsExtrap = 15; end
if ~exist('windowSize', 'var'), windowSize = 5; end
if ~exist('indEvery','var'), indEvery = 25; end
%-- For the toy2 data
if ~exist('rbfVariance', 'var'), rbfVariance = 0.1; end
if ~exist('outputNoise', 'var'), outputNoise = 0.0001; end
%--
if ~exist('initVardistIters', 'var'), initVardistIters = [400 300 0 0]; end
if ~exist('itNo', 'var'), itNo = [0 0 800 400 400]; end
if ~exist('propUncertIters', 'var'), propUncertIters = 50; end
if ~exist('resultsFile', 'var'), resultsFile = ['results/KstepAhead/' dataSetName '/results_Exp' num2str(experimentNo) '.txt']; end
if ~exist('resultsMatFile', 'var'), resultsMatFile = ['results/KstepAhead/' dataSetName '/results_Exp' num2str(experimentNo) '.mat']; end
if ~exist('diaryFile', 'var'), diaryFile = ''; end
if ~exist('dataSetName', 'var'), dataSetName = 'toy'; end % Other: mgdata, sunspot,
if ~exist('overWriteResults', 'var'), overWriteResults = false; end
if ~isempty(diaryFile)
    try
        delete(diaryFile)
    end
    diary(diaryFile);
    fprintf('\n\n\n');
    fprintf('############################################\n');
    fprintf('#   EXPERIMENT NO: %d                      #\n', experimentNo);
    fprintf('############################################\n\n');
end


vargplvm_init;

Ntr = windowsTr * windowSize;
Ninterp = windowsInterp * windowSize;
Nextrap = windowsExtrap * windowSize;
N = Ntr + Nextrap + windowSize;

switch dataSetName
    case 'mgdata'
        Yorig = vargplvmLoadData('mgdata');
    case 'sunspot'
        YY = sunspotData();
        Yorig = YY(:,2); clear('YY');
    case 'toy1'
        NN = 1000;
        rep = 20;
        t = linspace(0, rep*2*pi, NN)';
        kern = kernCreate(t, {'rbfperiodic', 'rbf'});%, 'lin'});
        kern.comp{1}.inverseWidth = 1;
        kern.comp{2}.variance = 0.1;
        kern.comp{2}.inverseWidth = 5;
        %kern.comp{3}.variance = 0.0000001;
        Kf = kernCompute(kern, t);
        Yorig = gsamp(zeros(NN, 1), Kf, 1)';
        Yorig = Yorig + sqrt(outputNoise).*randn(size(Yorig));
    case 'toy2'
        rep = 30;
        NN = 1000;
        t = linspace(0, rep*2*pi, NN)';
        kern = kernCreate(t, {'rbfperiodic', 'rbf', 'lin'});
        kern.comp{1}.inverseWidth = 1;
        kern.comp{2}.variance = rbfVariance;
        kern.comp{2}.inverseWidth = 0.5;
        kern.comp{3}.variance = 0.00015;
        Kf = kernCompute(kern, t);
        Yorig = gsamp(zeros(NN, 1), Kf, 1)'; % Change 1 to something else for multivariate!
        Yorig = Yorig + sqrt(outputNoise).*randn(size(Yorig));
        Yorig = Yorig(1:min(N,NN), :);
end





Yall= Yorig;
Yorig = Yorig(1:N,:);meanYorig = mean(Yorig(:));
Yorig = Yorig - mean(Yorig(:));
stdYall1 = std(Yorig);
Yorig = Yorig./std(Yorig);

Yall = Yall - meanYorig;
Yall = Yall./stdYall1;


% The real dataset is [y1,y2,.., y_K | y_{K+1}], ...
% The X indices are just for the plots...
[X, Y] = util_transformTimeSeriesToSeq(Yorig, windowSize);
%assert(sum(sum(abs(Yorig - util_transformSeqToTimeSeries(X, Y, windowSize)))) == 0);


numWindows = Ntr / windowSize;
perm = randperm(numWindows-2)+1; % Exclude 1st and last windows
interpWindows = perm(1:windowsInterp);
extrapWindows = numWindows + 1:numWindows+windowsExtrap;

% Remove some windows from tr. set for extrapolation
indInterp = [];
for i=sort(interpWindows)
    indInterp = [indInterp (i-1)*windowSize+1:(i-1)*windowSize+windowSize];
end

indExtrap = [];
for i=sort(extrapWindows)
    indExtrap = [indExtrap (i-1)*windowSize+1:(i-1)*windowSize+windowSize];
end

indTr = setdiff(1:size(X,1), [indInterp indExtrap]);


Xtr = X(indTr, :);
Ytr = Y(indTr, :);
Xinterp = X(indInterp, :);
Yinterp = Y(indInterp, :);
Xextrap = X(indExtrap, :);
Yextrap = Y(indExtrap, :);

t_tr = indTr';
t_interp = indInterp';
t_extrap = indExtrap';

%--------- Larger dataset
[X2,Y2]=util_transformTimeSeriesToSeq(Yall, windowSize);
numWindows = Ntr / windowSize;
extrapWindows = numWindows + 1:numWindows+windowsExtrap;
% Remove some windows from tr. set for extrapolation
indInterp = [];
for i=sort(interpWindows)
    indInterp = [indInterp (i-1)*windowSize+1:(i-1)*windowSize+windowSize];
end
indExtrap = [];
for i=sort(extrapWindows)
    indExtrap = [indExtrap (i-1)*windowSize+1:(i-1)*windowSize+windowSize];
end
YextrapMore = Y2(indExtrap(1):end,:);
XextrapMore = X2(indExtrap(1):end,:);
%--------------


%
% plot(t_tr, Yorig(t_tr, 1), 'x', 'MarkerSize', 6, 'LineWidth', 2); hold on
% plot(t_interp, Yorig(t_interp, 1), 'Ok', 'MarkerSize', 6, 'LineWidth', 2);
% plot(t_extrap, Yorig(t_extrap, 1), 'sr', 'MarkerSize', 6, 'LineWidth', 2);
% plot(1:size(Yorig,1), Yorig(:,1), '-g');
% legend({'Training', 'inter', 'extrap'})
% hold off
%%





%% BGPLVM Uncertain with NOT fixed inducing and fixed X.
fileName = vargplvmWriteResult([], 'vargplvm', dataSetName, experimentNo, '', './matFiles/KstepAhead/');
options_vargplvm = vargplvmOptions('dtcvar');
options_vargplvm.kern = globalOpt.mappingKern;
options_vargplvm.numActive = size(Ytr,1);
options_vargplvm.optimiser = 'scg2';
options_vargplvm.initSNR = 100;

options = options_vargplvm;
options.numActive = min(options_vargplvm.numActive, 30);
options.fixInducing = false;
options.learnInducing = true;
options.initX = Xtr;
options.fixVardist = true;
if exist([fileName '.mat'], 'file') && ~overWriteResults
    load(fileName);
    fprintf('# Found and loaded file %s on disk!\n', fileName);
    model = vargplvmRestorePrunedModel(modelPruned, Ytr);
else
    model = vargplvmCreate(size(Xtr,2), size(Ytr,2), Ytr, options);
    model = vargplvmModelInit(model, globalOpt);
    model.vardist.covars = 1e-10*model.vardist.covars;  % Very small variance for tr. data
    % Do not learn training X's (should be done better, ie not calculate gradients of X's at all)
    model.windowSize = windowSize;
    if exist('parallel','var') && parallel && (isfield(options,'learnInducing') && options.learnInducing),model.vardist.parallel = true;  end
    model = vargplvmOptimiseModel(model, 0, 0, {globalOpt.initVardistIters, globalOpt.itNo});
end


if ~exist([fileName '.mat'], 'file') || overWriteResults
    modelPruned = vargplvmPruneModel(model);
    save([fileName '.mat'], 'modelPruned');
end


if ~isempty(diaryFile)
    diary off;
end


%% PREDICTIONS (different approaches)


%------------- vargplvm without adding Ind pts ------------
opt.k=1;                 opt.augmentModel=true; opt.augmentInd=false;
opt.trCovarsMult = 1;    opt.varInPred=true;    opt.reOptimiseIters=0;
% Uncomment for quicker if you dont want the large extrapolation
%[Z3, varZ3] = vargplvmPredictUncertainK(model, Xextrap(1,:), size(Xextrap, 1),opt);
if ~exist('Z3more', 'var') || ~exist('varZ3more','var')
    [Z3more, varZ3more,m3] = vargplvmPredictUncertainK(model, XextrapMore(1,:), size(XextrapMore, 1),opt);
end
Z3 = Z3more(1:size(Yextrap,1),:); varZ3 = varZ3more(1:size(Yextrap,1),:);
if ~exist('Z3itmore', 'var') || ~exist('varZ3itmore','var')
    opt.k_it = 100;   opt.ind_it = 15; opt.augmentModel_it = 0;
    [Z3itmore, varZ3itmore]=vargplvmPredictUncertainWrap(model, XextrapMore, YextrapMore, opt);
end
Z3it = Z3itmore(1:size(Yextrap,1),:); varZ3it = varZ3itmore(1:size(Yextrap,1),:);


%-------------- Add noise ONLY in inputs and NOT in prediction.
opt.k=1;                 opt.augmentModel=true;  opt.augmentInd=false;
opt.trCovarsMult = 1;    opt.varInPred=false;    opt.reOptimiseIters=0;
% Uncomment for quicker if you dont want the large extrapolation
%[Z33, varZ33] = vargplvmPredictUncertainK(model, Xextrap(1,:), size(Xextrap, 1),opt);
if ~exist('Z33more', 'var') || ~exist('varZ33more','var')
    [Z33more, varZ33more,m33] = vargplvmPredictUncertainK(model, XextrapMore(1,:), size(XextrapMore, 1),opt);
end
Z33 = Z33more(1:size(Yextrap,1),:); varZ33 = varZ33more(1:size(Yextrap,1),:);


%------------- vargplvm with adding Ind pts at test time ------------
% Old results (better for predictions, worse for uncertainty - also, slower) - equiv. to Zex_fixX2
opt.k=1;                 opt.augmentModel=true; opt.augmentInd=indEvery;
opt.trCovarsMult = 1;    opt.varInPred=true;    opt.reOptimiseIters=0;
%[Z3ind, varZ3ind] = vargplvmPredictUncertainK(model, Xextrap(1,:), size(Xextrap, 1),opt);
if ~exist('Z3moreInd', 'var') || ~exist('varZ3moreInd','var')
    [Z3moreInd, varZ3moreInd,m3Ind] = vargplvmPredictUncertainK(model, XextrapMore(1,:), size(XextrapMore, 1),opt);
end
Z3ind = Z3moreInd(1:size(Yextrap,1),:); varZ3ind = varZ3moreInd(1:size(Yextrap,1),:);


savePlots = false;
saveDir = 'autoreg/';
axFs = 18;
titleFs = 0;
legFs = 18;



%% Evaluation
if ~exist('evaluation','var'), evaluation = false; end

if evaluation
    % Errors
    MAE_Z3 = mean(abs(Z3more(:) - YextrapMore(:)));
    MSE_Z3 = mean((Z3more(:) - YextrapMore(:)).^2);

    fprintf('######## ERRORS #############\n\n');
    fprintf('# MAE: %.3f \n', MAE_Z3);
    fprintf('# MSE: %.3f \n', MSE_Z3);
    
    %
    %-- Plots
    ymax = max([max(Z3more) max(YextrapMore)]);
    ymin = min([min(Z3more) min(YextrapMore)]);
    space = (ymax-ymin)./18;
    
    t_extrap_more = [t_tr(end)+1:t_tr(end)+size(YextrapMore,1)]';
    close all;
    figs = {};
    figs{end+1} = figure;
    
    lgnd = {};
    names={'true','ours'}
    symb={'k-+','bs-'}
    % Ours
    plot(t_extrap_more, YextrapMore, symb{1});  lgnd{end+1}=names{1}; hold on; grid on;
    plot(t_extrap_more, Z3more, symb{2});       lgnd{end+1}=names{2};
    set(gca, 'XtickLabel',[]);
    set(gca, 'FontSize', axFs);
    xlim([t_extrap_more(1)-20 t_extrap_more(end)+10]);
    ylim([ymin - space ymax+space])
    
    % Variances
    figs{end+1} = figure; lgnd = {};
    plot(t_extrap_more, 2*sqrt(varZ3more), symb{2});       lgnd{end+1}=names{2}; hold on; grid on
    set(gca, 'FontSize', axFs);
    xlim([t_extrap_more(1)-20 t_extrap_more(end)+10]);
    ylim([0 0.13])
    
end


