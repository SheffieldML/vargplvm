% DEMCMU35VARGPLVMPLOTSSCALED Load results for the dyn. GPLVM on CMU35 data and produce plots
% in the scaled space (compatible with FGPLVM plots)
% DESC Load results for the dyn. GPLVM on CMU35 data and produce plots in
% the scaled space.
% COPYRIGHT :  Andreas C. Damianou, Michalis K. Titsias, 2011
%
% SEEALSO : demCmu35VargplvmLoadChannels, demCmu35gplvmVargplvm3.m,
%           demCmu35VargplvmAnimate.m
% VARGPLVM



randn('seed', 1e5);
rand('seed', 1e5);


dataSetName = 'cmu35gplvm';
if ~exist('experimentNo') experimentNo = 33; end
% Options for the following: 'Legs', 'Body'
if ~exist('predictPart') predictPart = 'Legs'; end
% Sampling:
if ~exist('doSampling') doSampling = 0; end
if ~exist('showSkel') showSkel = 1; end
if ~exist('displayPlots')  displayPlots = 1; end
% -1 is a flag meaning that plotRange == missingInd
if ~exist('plotRange') plotRange = -1; end

%fprintf(1,'# Preparing dataset and loading results...\n');

% load data
[Y, lbls, Ytest, lblstest] = lvmLoadData(dataSetName);
dataSetName = 'cmu35gplvm';

capName = dataSetName;
capName(1) = upper(capName(1));
fileName=['dem' capName 'Vargplvm' num2str(experimentNo)];

% load data
[Y, lbls, Ytest, lblstest] = lvmLoadData(dataSetName);
origBias = mean(Y);
origScale = 1./sqrt(var(Y));

legInd = [8:14];
bodyInd = [21:50];

startInd = 63;
dt=0.05;
timeStampsTest = ([0:size(Ytest,1)-1].*dt)';

load(fileName)
model.dynamics.t_star = timeStampsTest;

% Choose one of the following lines as appropriate to obtain the results for the leg or body reconstruction
% respectively.
if strcmp(predictPart,'Legs')
    load(['demCmu35Vargplvm' num2str(experimentNo) 'PredLegs.mat']); missingInd = legInd;
elseif strcmp(predictPart, 'Body')
    load(['demCmu35Vargplvm' num2str(experimentNo) 'PredBody.mat']) ;missingInd = bodyInd;
end

%------
Ygplvm=Y;
YtestGplvm=Ytest;
YtestGplvmOrig = Ytest;
Y = Y - repmat(origBias, size(Y, 1), 1);
Ytest = Ytest - repmat(origBias, size(Ytest, 1), 1);
Y = Y.*repmat(origScale, size(Y, 1), 1);
Ytest = Ytest.*repmat(origScale, size(Ytest, 1), 1);

origYtest = Ytest;
YtrueTest = Ytest;

temp = ones(1, size(Ytest, 2));
temp(missingInd) = 0;
presentInd = find(temp);

YtrainNn = Y(:, presentInd);
YtestNn = origYtest(startInd:end, presentInd);
%---------




YtestGplvm(startInd:end, missingInd) = NaN;
indexMissingData = startInd:size(YtestGplvm);
% In case only barmu and lambda are loaded from the results then the next step is necessary to find
% the predictive dencity Ypred.
[x, varx, modelUpdated] = vargplvmDynamicsUpdateModelTestVar(model, barmu, lambda, YtestGplvm);
Testmeans = x(indexMissingData, :);
Testcovars = varx(indexMissingData, :);
[mu, sigma] = vargplvmPosteriorMeanVar(modelUpdated, Testmeans, Testcovars);
Ypred = mu;


%%
% NN
% First nearest neighbour
dists = dist2(YtrainNn./(repmat(origScale(presentInd), size(YtrainNn, 1), 1)), ...
              YtestNn./(repmat(origScale(presentInd), size(YtestNn, 1), 1)));
[void, bestIndNn] = min(dists);
lenVal = size(Ytest, 1);
err = (YtrueTest(startInd:end, missingInd) - Y(bestIndNn, missingInd))...
         ./repmat(origScale(1, missingInd), size(YtrueTest, 1)-startInd+1, 1);
err = err*180/pi;
errStruct.angleErrorNnScaled = sqrt(mean(mean((err.*err))));
load cmu35TaylorScaleBias
err = err.*repmat(scale(1, missingInd+9), size(YtrueTest, 1)-startInd+1, 1);
errStruct.taylorErrorNnScaled = sum(sum(err.*err))/length(missingInd);

% Second nearest neighbour
dists = dist2(YtrainNn, ...
              YtestNn);
[void, bestIndNn] = min(dists);
lenVal = size(Ytest, 1);
err = (YtrueTest(startInd:end, missingInd) - Y(bestIndNn, missingInd))...
         ./repmat(origScale(1, missingInd), size(YtrueTest, 1)-startInd+1, 1);
err = err*180/pi;
errStruct.angleErrorNn = sqrt(mean(mean((err.*err))));
load cmu35TaylorScaleBias
err = err.*repmat(scale(1, missingInd+9), size(YtrueTest, 1)-startInd+1, 1);
errStruct.taylorErrorNn = sum(sum(err.*err))/length(missingInd);


%%


YtrueTest=YtestGplvmOrig;
err = YtrueTest(startInd:end, missingInd) - Ypred(:, missingInd);

% Convert to degrees.
err = err*180/pi;
% Compute average of mean square error.
errStruct.angleErrorGplvm = sqrt(mean(mean((err.*err))));
load cmu35TaylorScaleBias
err = err.*repmat(scale(1, missingInd+9), size(YtrueTest, 1)-startInd+1, 1);
errStruct.taylorErrorGplvm = sum(sum(err.*err))/length(missingInd);


 % Normalise predictions as in Fgplvm  
 Ypred = Ypred - repmat(origBias, size(Ypred, 1), 1);
 Ypred = Ypred.*repmat(origScale, size(Ypred, 1), 1);


if plotRange == -1
    plotRange = missingInd;
end
colordef white
for plotNo = plotRange
    figNo = plotNo - min(plotRange) + 1;
    figure(figNo)
    clf
    lin = plot(1:size(origYtest, 1), origYtest(:, plotNo), '-');
    hold on
    lin = [lin; plot(1:size(origYtest, 1), [origYtest(1:startInd-1, plotNo); Y(bestIndNn, plotNo)], ':')];
    lin = [lin; plot(1:size(origYtest, 1), [origYtest(1:startInd-1, plotNo); Ypred(:, plotNo)], '--')];
    ax = gca;
    xlabel('Frame')
    ylabel('Normalized joint angle');
    set(ax, 'fontname', 'arial');
    set(ax, 'fontsize', 20);
    set(lin, 'lineWidth', 2);
    
    % make smaller for PNG plot.
    pos = get(gcf, 'paperposition');
    origpos = pos;
    pos(3) = pos(3)/2;
    pos(4) = pos(4)/2;
    set(gcf, 'paperposition', pos);
    fontsize = get(gca, 'fontsize');
    set(gca, 'fontsize', fontsize/2);
    lineWidth = get(gca, 'lineWidth');
    set(gca, 'lineWidth', lineWidth*2);
    %print('-dpng', ['../html/' fileName])
    %plot2svg(['../html/' fileName '.svg'])
    set(gcf, 'paperposition', origpos);
end

