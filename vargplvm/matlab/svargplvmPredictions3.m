% SVARGPLVMPREDICTIONS3 This script implements the predictions for a classification task
%
% SEEALSO : demClassification, demClassification3, demClassificationGeneral, svargplvmPredictions2, svargplvmPredictions
%
% COPYRIGHT : Andreas C. Damianou, 2011

% VARGPLVM

%%%--- Script: svargplvmPredictions2 %%%%-
%%% This script is like svargplvmPredictions with a few additions (it
%%% should replace svargplvmPRedictions)
if ~exist('testOnTraining')
    testOnTraining=0;
end
if ~exist('globalOpt') || ~isfield(globalOpt, 'predictionsOptionX')
    predictionsOptionX = 1;
end
if ~exist('globalOpt') || ~isfield(globalOpt, 'predictionsGenerationOption')
    predictionsGenerationOption = 1;
end
if ~exist('testDisplayIters')
    testDisplayIters=0; %%%
end
if ~exist('predictionVariance')
    predictionVariance = false;
end

numberTestPoints = size(Yts{1},1); % 10;
if testOnTraining
   % perm = randperm(model.N);
   % testInd = perm(1:size(Yall{1},1));
   testInd = 1:size(Yall{1},1);
else
    perm = randperm(size(Yts{obsMod},1));
    %testInd = perm(1:numberTestPoints);
    if ~exist('testInd')
        testInd = 1:size(Yts{1},1); %%%%  
    end
end


ZpredMuAll = [];
XpredAll = zeros(length(testInd), model.q);  %%%
varZpredAll = zeros(length(testInd), model.q);  %%%
pb = myProgressBar(length(testInd));
for i=1:length(testInd)
    pb = myProgressBar(pb,i);
    curInd = testInd(i);
   % fprintf('# Testing indice number %d ', curInd);
    if testOnTraining
        y_star = model.comp{obsMod}.y(curInd,:);
        x_star = model.comp{obsMod}.vardist.means(curInd,:);
        varx_star = model.comp{obsMod}.vardist.covars(curInd,:);
    else
        y_star = Yts{obsMod}(curInd,:);
        z_star = Yts{infMod}(curInd,:);
        dst = dist2(y_star, model.comp{obsMod}.y);
        [mind, mini(i)] = min(dst);
        
        model.comp{obsMod}.mini = mini(i);
        Init(i,:) = model.vardist.means(mini(i),:);
        vardistx = vardistCreate(model.comp{obsMod}.vardist.means(mini(i),:), model.q, 'gaussian');
        vardistx.covars = model.comp{obsMod}.vardist.covars(mini(i),:);
        model.comp{obsMod}.vardistx = vardistx;
        iters = globalOpt.reconstrIters;
        % Find p(X_* | Y_*) which is approximated by q(X_*)
        [x_star, varx_star, modelUpdated] = vargplvmOptimisePoint(model.comp{obsMod}, vardistx, y_star, testDisplayIters, iters);%%%
    end
    numberOfNN = 1;
    % Now we selected a datapoint X_* by taking into account only the
    % private dimensions for Y. Now, based on the shared dimensions of
    % that, we select the closest (in a NN manner) X from the training data.
   % fprintf('# Finding the %d NN of X_* with the training X based only on the shared dims.\n', numberOfNN);
    [ind, distInd] = nn_class(model.X(:,sharedDims), x_star(:,sharedDims), numberOfNN, 'euclidean');
    
    ZpredMu = zeros(length(ind), size(model.comp{infMod}.y,2));
    varZstar = zeros(size(ZpredMu));
    ZpredSigma = zeros(length(ind), size(model.comp{infMod}.y,2));
      
    
    % Find p(y_*|x_*) for every x_* found from the NN
   % fprintf('# Predicting from the NN of X_* ');
    for k=1:numberOfNN
        x_cur = svargplvmPredictX(x_star, model.X, ind(k), sharedDims, predictionsOptionX);
        %[ZpredMu(k,:), ZpredSigma(k,:)] =
        %vargplvmPosteriorMeanVar(model.comp{infMod}, model.X(ind(k),:));
        if predictionsGenerationOption
            % NEW: varx_star
            if ~predictionVariance
                ZpredMu(k,:) = vargplvmPosteriorMeanVar(model.comp{infMod},x_cur, varx_star);
            else
                [ZpredMu(k,:) varZstar(k,:)] = vargplvmPosteriorMeanVar(model.comp{infMod},x_cur, varx_star);
            end
        else
            [ind2, distInd] = nn_class(model.X, x_cur, 1,'euclidean');
            ZpredMu(k,:) = model.comp{infMod}.y(ind2,:);
        end
    end
    ZpredMuAll = [ZpredMuAll; ZpredMu(1,:)];
    if predictionVariance
        varZpredAll = [varZpredAll; varZstar(1,:)];
    end
    XpredAll(i,:) = x_cur;
   % fprintf('\n\n');
end
fprintf(1, '\n');

ZpredMuAllOrig = ZpredMuAll;
defThresh = 0.5;
if min(min(Yts{2} == -1))
    defThresh = 0;
end
ZpredMuAll(ZpredMuAll >= defThresh) = 1;
ZpredMuAll(ZpredMuAll < defThresh) = 0;


