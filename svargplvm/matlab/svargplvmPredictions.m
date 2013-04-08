% SHEFFIELDMLPREDICTIONS This script implements the predictions for a classification task
%
% SEEALSO : demClassification, demClassification3, demClassificationGeneral, svargplvmPredictions2, svargplvmPredictions3
%
% COPYRIGHT : Andreas C. Damianou, 2011

% SHEFFIELDML


%%%--- Script: svargplvmPredictions %%%%-

numberTestPoints = size(Yts{1},1); % 10;
if testOnTraining
    perm = randperm(model.N);
    testInd = perm(1:numberTestPoints);
else
    perm = randperm(size(Yts{obsMod},1));
    %testInd = perm(1:numberTestPoints);
    if ~exist('testInd')
        testInd = 1:size(Yts{1},1); %%%%  
    end
end


ZpredMuAll = [];
XpredAll = zeros(length(testInd), model.q);  %%%
varXpredAll = []; % zeros(length(testInd), model.q);  %%%
pb = myProgressBar(length(testInd), min(length(testInd),20));
for i=1:length(testInd)
    pb = myProgressBar(pb,i);
    curInd = testInd(i);
   % fprintf('# Testing indice number %d ', curInd);
    if testOnTraining
     %   fprintf('taken from the training set\n');
        y_star = model.comp{obsMod}.y(curInd,:);
        x_star = model.comp{obsMod}.vardist.means(curInd,:);
        varx_star = model.comp{obsMod}.vardist.covars(curInd,:);
    else
      %  fprintf('taken from the test set\n');
        y_star = Yts{obsMod}(curInd,:);
        z_star = Yts{infMod}(curInd,:);
        dst = dist2(y_star, model.comp{obsMod}.y);
        [mind, mini(i)] = min(dst);
        
        model.comp{obsMod}.mini = mini(i);
        Init(i,:) = model.vardist.means(mini(i),:);
        vardistx = vardistCreate(model.comp{obsMod}.vardist.means(mini(i),:), model.q, 'gaussian');
        vardistx.covars = model.comp{obsMod}.vardist.covars(mini(i),:);
        model.comp{obsMod}.vardistx = vardistx;
        display=0; %%%
        iters = globalOpt.reconstrIters;
        % Find p(X_* | Y_*) which is approximated by q(X_*)
        [x_star, varx_star, modelUpdated] = vargplvmOptimisePoint(model.comp{obsMod}, vardistx, y_star, display, iters);%%%
    end
    numberOfNN = 1;
    % Now we selected a datapoint X_* by taking into account only the
    % private dimensions for Y. Now, based on the shared dimensions of
    % that, we select the closest (in a NN manner) X from the training data.
   % fprintf('# Finding the %d NN of X_* with the training X based only on the shared dims.\n', numberOfNN);
    [ind, distInd] = nn_class(model.X(:,sharedDims), x_star(:,sharedDims), numberOfNN, 'euclidean');
    
    ZpredMu = zeros(length(ind), size(model.comp{infMod}.y,2));
    ZpredSigma = zeros(length(ind), size(model.comp{infMod}.y,2));
      
    
    % Find p(y_*|x_*) for every x_* found from the NN
   % fprintf('# Predicting from the NN of X_* ');
    for k=1:numberOfNN
        %      fprintf('.');
        x_cur = model.X(ind(k),:);
        %x_cur(sharedDims) = x_star(sharedDims); %%% OPTIONAL!!!
        %[ZpredMu(k,:), ZpredSigma(k,:)] = vargplvmPosteriorMeanVar(model.comp{infMod}, model.X(ind(k),:));

        %ZpredMu(k,:) = vargplvmPosteriorMeanVar(model.comp{infMod}, x_cur);

        [ind2, distInd] = nn_class(model.X, x_cur, 1,'euclidean'); %%-- NEW
        ZpredMu(k,:) = model.comp{infMod}.y(ind2,:); %%-- NEW


    end
    ZpredMuAll = [ZpredMuAll; ZpredMu(1,:)];
    varXpredAll = [varXpredAll; varx_star];
    XpredAll(i,:) = x_cur;
   % fprintf('\n\n');
end
fprintf(1, '\n');
