% DEMCMU35VARGPLVMANIMATE Load results for the dyn. GPLVM on CMU35 data and
% produce animation and plots. The plots (and errors) are in the data
% space. For plots / errors  in the scaled space (for comparison with FGPLVM) 
% see demCmu35VargplvmPlotsScaled.m
% DESC Load results for the dyn. GPLVM on CMU35 data and produce animation
% and plots.
% COPYRIGHT :  Andreas C. Damianou, Michalis K. Titsias, 2011
%
% SEEALSO : demCmu35VargplvmLoadChannels, demCmu35gplvmVargplvm3.m,
% demCmu35VargplvmPlotsScaled.m
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

YtestOrig = Ytest;
YtestGplvm = Ytest;
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
origYtest = Ytest;
err = origYtest(startInd:end, missingInd) - Ypred(:, missingInd);
err = err*180/pi;
angleErrorGplvm = sqrt(mean(mean((err.*err))))
load cmu35TaylorScaleBias
err = err.*repmat(scale(1, missingInd+9), size(origYtest, 1)-startInd+1, 1);
taylorErrorGplvm = sum(sum(err.*err))/length(missingInd)


%%
% NN
temp = ones(1, size(Ytest, 2));
temp(missingInd) = 0;
presentInd = find(temp);
YtrainNn = Y(:, presentInd);
YtestNn = origYtest(startInd:end, presentInd);
dists = dist2(YtrainNn, YtestNn);
[void, bestIndNn] = min(dists);

%%



if displayPlots
    colordef white
    if plotRange == -1 % flag
        plotRange = missingInd;
    end
    for plotNo = plotRange
        figNo = plotNo - min(plotRange) + 1;
        figure(figNo)
        clf
        lin = plot(1:size(origYtest, 1), origYtest(:, plotNo), '-')
        hold on
        lin = [lin; plot(1:size(origYtest, 1), [origYtest(1:startInd-1, plotNo); Y(bestIndNn, plotNo)], ':')];
        lin = [lin; plot(1:size(origYtest, 1), [origYtest(1:startInd-1, plotNo); Ypred(:, plotNo)], '--')];
        ax = gca;
        xlabel('Frame')
        ylabel('Normalized joint angle')
        set(ax, 'fontname', 'arial');
        set(ax, 'fontsize', 20);
        set(lin, 'lineWidth', 2);
        pos = get(gcf, 'paperposition')
        origpos = pos;
        pos(3) = pos(3)/2;
        pos(4) = pos(4)/2;
        set(gcf, 'paperposition', pos);
        fontsize = get(gca, 'fontsize');
        set(gca, 'fontsize', fontsize/2);
        lineWidth = get(gca, 'lineWidth');
        set(gca, 'lineWidth', lineWidth*2);
        set(gcf, 'paperposition', origpos);
    end
end


%%

if showSkel
    skel = acclaimReadSkel('35.asf');
    [tmpchan, skel] = acclaimLoadChannels('35_01.amc', skel);
    channels1 = demCmu35VargplvmLoadChannels(Ytest,skel);
    channels2 = demCmu35VargplvmLoadChannels(Ypred,skel);
    channels3 = demCmu35VargplvmLoadChannels(Y(bestIndNn, :),skel);
    %skelPlayData(skel, channels, 1/25);

    startInd = 63;
    skelPlayData2(skel, 1/15, channels1(startInd:end,:), channels2, channels3,{'Ytest','YpredGPLVM','YpredNN'});
end

%%



if doSampling
 %   timeStampsTestOrig = timeStampsTest;
 %   timeStampsTest = (timeStampsTest(end):dt:(timeStampsTest(end)+100*dt))';
    model.dynamics.t_star = timeStampsTest;
    seqTest = 1; % WALK
    %seqTest = 17; % RUN
    
    if seqTest ~= 1
        datFrom = model.dynamics.seq(seqTest-1)+1;
    else
        datFrom = 1;
    end
    datEnd=model.dynamics.seq(seqTest);
    
    model.dynamics.t_star = model.dynamics.t_star(1:(datEnd+1-datFrom));%%%
    %     clear model
    %     load demCmu35VargplvmOne
    %     %load demCmu35VargplvmOneRun
    %     ySamp, xSamp] = vargpTimeDynamicsSample(model,1);
    %     channels4 = demCmu35VargplvmLoadChannels(ySamp,skel);
    %     skelPlayData(skel, channels4, 1/25);
    
    
    [ySamp, xSamp] = vargpTimeDynamicsSampleSeq(model,1,seqTest);
    channels4 = demCmu35VargplvmLoadChannels(ySamp,skel);
    skelPlayData(skel, channels4, 1/10);
    
    %     [Testmeans2 Testcovars2] = vargplvmPredictSeq(model.dynamics, timeStampsTest,1);
    %     Varmu2 = vargplvmPosteriorMeanVar(model, Testmeans2, Testcovars2);
    %     channels5 = demCmu35VargplvmLoadChannels(Varmu2,skel);
    %     skelPlayData(skel, channels5, 1/25);
    
end
%%
%lvmVisualise(model, [], 'skelVisualise', 'skelModify', channels, skel);
%lvmVisualise(model, [], 'TEMPskelVisualise', 'skelModify', channels1,
%skel);