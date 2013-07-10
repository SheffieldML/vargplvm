function [errStruct] = vargplvmTaylorAngleErrors(model, Y, Ytest, startInd, origBias, origScale,...
    missingInd, name, expNo)


% VARGPLVMTAYLORANGLEERRORS Helper function for computing angle errors for CMU 35 data using
%	GPLVM with dynamics.
%
%	Description:
%	[errStruct] = vargplvmTaylorAngleErrors(model, Y, Ytest, startInd, origBias, origScale,...
%                                              missingInd, name, expNo);
% 	Based on fgplvmTaylorAngleErrors.m for FGPLVM.
% COPYRIGHT: Neil Lawrence, Andreas Damianou, Michalis Titsias 2010-2011
% VARGPLVM

%%%NEW
Ygplvm=Y;
YtestGplvm=Ytest;
YtestGplvmOrig = Ytest;

Y = Y - repmat(origBias, size(Y, 1), 1);
Ytest = Ytest - repmat(origBias, size(Ytest, 1), 1);
Y = Y.*repmat(origScale, size(Y, 1), 1);
Ytest = Ytest.*repmat(origScale, size(Ytest, 1), 1);
%%%


origYtest = Ytest;
YtrueTest = Ytest;

capName = name;
capName(1) = upper(capName(1));

temp = ones(1, size(Ytest, 2));
temp(missingInd) = 0;
presentInd = find(temp);

YtrainNn = Y(:, presentInd);
YtestNn = origYtest(startInd:end, presentInd);



disp('# NN prediction...');

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


Ytest(startInd:end, missingInd) = NaN;


%%%%%%%%%%%%%  VAR  GPLVM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
YtestGplvm(startInd:end, missingInd) = NaN;

disp('# Initializing latent points...');
indexPresent = setdiff(1:model.d, missingInd);

%iters=2000;% default: 1000
iters = model.reconstrIters;
fprintf(1,'# Optimising for %d iterations...\n',iters)

if isfield(model, 'dynamics') & ~isempty(model.dynamics)
    for i=1:size(YtestGplvm,1)
        % initialize the latent points using the nearest neighbour from the training data
        dst = dist2(YtestGplvm(i,indexPresent), Ygplvm(:,indexPresent));
        [mind, mini(i)] = min(dst);
    end
    
    
    
    
    indexMissingData = startInd:size(YtestGplvm);
    
    vardistx = vardistCreate(model.dynamics.vardist.means(mini,:), model.q, 'gaussian');
    vardistx.covars = 0.2*ones(size(vardistx.covars));
    model.vardistx = vardistx;
    
    % Do also reconstruction in test data
    [x, varx] = vargplvmOptimiseSeqDyn(model, vardistx, YtestGplvm, 1, iters);
    % keep the optimized variational parameters
    barmu = x;
    lambda = varx;
    
    % Get the variational means and variacnes for the new test sequcen and
    % update the model to be prepared for prediction
    [x, varx, modelUpdated] = vargplvmDynamicsUpdateModelTestVar(model, barmu, lambda, YtestGplvm);
    
    % Latent variables corresponding to the data  with missing dimensions
    Testmeans = x(indexMissingData, :);
    Testcovars = varx(indexMissingData, :);
    
    %  YtsOrMs = YtsOriginal(indexMissingData,:);
    disp('# Making the predictions for Y...')
    % Reconstruct all dimensions (missing and observed) for the test data
    [mu, sigma] = vargplvmPosteriorMeanVar(modelUpdated, Testmeans, Testcovars);
    Ypred = mu;
else
    %... STatic...
    % Static case doesnt capture data dependencies. Hence, we reconstruct
    % only the missing inputs:
    YtestGplvm = YtestGplvm(startInd:end, :);

    
   for i=1:size(YtestGplvm,1)
        % initialize the latent points using the nearest neighbour from the training data
        dst = dist2(YtestGplvm(i,indexPresent), Ygplvm(:,indexPresent));
        [mind, mini(i)] = min(dst);
    end

   
    vardistx = vardistCreate(model.vardist.means(mini,:), model.q, 'gaussian');
    vardistx.covars = 0.2*ones(size(vardistx.covars));
    model.vardistx = vardistx;
    % mu=[];
    % sigma=[];
    % for i=1:size(YtestGplvm,1)
    % optimize mean and vars of the latent point
    model.vardistx = vardistx;
    
    [x, varx] = vargplvmOptimisePoint(model, vardistx, YtestGplvm, 1, iters); %
    
    % reconstruct the missing outputs
    [mu, sigma] = vargplvmPosteriorMeanVar(model, x, varx);
    Ypred = mu;
    %    mu(i,:) = mu_i;
    %    sigma(i,:) = sigma_i;
    %end
    %save 'temp_temp.mat' %%%%
end

%{
[Xpred, varx] = vargplvmOptimiseSeqDyn(model, vardistx, YtestGplvm, 1, iters); % default: 1000
save('TEMPtaylorVargplvm2.mat','model','vardistx');
disp('# Getting the original variational parameters...');
barmu = vardistx.means;
lambda = vardistx.covars;
modelTest = vargplvmDynamicsUpdateModelTestVar(model, barmu, lambda);
Nstar = size(YtestGplvm,1);
overL=0;
Xpred = modelTest.vardist.means(end-(Nstar-1+overL):end,:);
varx= modelTest.vardist.covars(end-(Nstar-1+overL):end,:);
%}


%%%%%%%%%%%%%
% fName=['TEMPTaylorVargplvm' num2str(iters) 'it' num2str(expNo) '.mat'];
% fprintf(1,'# Saving the workspace in file %s...\n',fName);
% save(fName);

%%%%%%%%%


%%% Normalise predictions as in Fgplvm
% Ypred = Ypred - repmat(origBias, size(Ypred, 1), 1);
% Ypred = Ypred.*repmat(origScale, size(Ypred, 1), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




disp('# Final calculations...')

YtrueTest=YtestGplvmOrig;

% % De-normalise the predictions.
% err = (YtrueTest(startInd:end, missingInd) - Ypred(:, missingInd))./repmat(origScale(1, missingInd), size(YtrueTest, 1)-startInd+1, 1);

err = YtrueTest(startInd:end, missingInd) - Ypred(:, missingInd);

% Convert to degrees.
err = err*180/pi;
% Compute average of mean square error.
errStruct.angleErrorGplvm = sqrt(mean(mean((err.*err))));
load cmu35TaylorScaleBias
err = err.*repmat(scale(1, missingInd+9), size(YtrueTest, 1)-startInd+1, 1);
errStruct.taylorErrorGplvm = sum(sum(err.*err))/length(missingInd);


%%% Normalise predictions as in Fgplvm  %NEW!!!!!!!!!!
Ypred = Ypred - repmat(origBias, size(Ypred, 1), 1);
Ypred = Ypred.*repmat(origScale, size(Ypred, 1), 1);
%%%%

plotRange = missingInd;
colordef white
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
    fileName = ['dem' capName 'Reconstruct' num2str(expNo) '_' num2str(plotNo)];
    %  print('-depsc', ['../tex/diagrams/' fileName])
    %  print('-deps', ['../tex/diagrams/' fileName 'NoColour'])
    %plot2svg(['../html/' fileName '.svg'])
    
    % make smaller for PNG plot.
    pos = get(gcf, 'paperposition')
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


%%% Denormalise to save:  %NEW!!!!!!!!!!
Ypred = Ypred + repmat(origBias, size(Ypred, 1), 1);
Ypred = Ypred./repmat(origScale, size(Ypred, 1), 1);
%%%%

%%%
%prunedModelUpdated = vargplvmPruneModel(modelUpdated,1);
fileToSave = 'demCmu35Vargplvm';
if strcmp(name(end-3:end),'Body')
    n = 'PredBody.mat';
else
    n='PredLegs.mat';
end

if isfield(model, 'dynamics') & ~isempty(model.dynamics)
	prunedModel = vargplvmPruneModel(model);
	prunedModelUpdated = vargplvmPruneModel(modelUpdated);
	save([fileToSave num2str(expNo) n], 'prunedModel','prunedModelUpdated','Ypred','errStruct');
else
	save([fileToSave num2str(expNo) n], 'model','Ypred','errStruct');
end
fprintf(1,'# Saved %s\n',[fileToSave num2str(expNo) n]);
%%%
