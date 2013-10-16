% SVARGPLVMPREDICTIONS This script implements the predictions for MRD
%
% SEEALSO : demClassification, demClassification3, demClassificationGeneral, svargplvmPredictions2, svargplvmPredictions3
%
% COPYRIGHT : Andreas C. Damianou, 2011

% VARGPLVM

% There are two kinds of predictions: on training and on test data.
% * Predict on the training data: given a test points y* which is a copy of
% a training point y, we want to find the K points z*_1, z*_2, ..., z*_K in
% the other modality that best match y in the latent space, ie if output y
% corresponds to a latent point x_y, we want to find the K latent points
% x_z1, x_z2, ... that are most similar to x_y (in the shared dimensions
% only!!) from which we generate the corresponding z*s (you select any K
% you want). This is how we
% solve the corresponding problem (see paper, and Yale faces demo). Of
% course, given a y* from the training points (e.g. Ytr{1}(i,:), the most likely z* is going
% to be the z* of the training data on the other modality (Ytr{2}(i,:)),
% but the important thing is to see which other points are returned.
%
% * Predict in the test data: as above, but the test point y* is totally
% unobserved (new). Then, we have to learn the corresponding x*_y by
% optimising a variational distribution, and generate K (for some K) points
% in the other modality as explained in the paper. See also the
% classification demo and the code in the current script.

%%%--- Script: svargplvmPredictions %%%%-

% This is the number that is referred to as K in the top of this
% script, ie how many points in the other modality we want to recover /
% generate.
if ~exist('numberOfNN', 'var'), numberOfNN = 1; end
if ~exist('saveAllNN', 'var'), saveAllNN = 0; end
if ~exist('infMethod','var'), infMethod = 1; end
if ~exist('displayTestOpt', 'var'), displayTestOpt=false; end
% TODO:
% Here replace model.comp{i}.m with model.comp{i}.mOrig if DgtN is
% active...


if testOnTraining
    if ~exist('numberTestPoints', 'var')
        numberTestPoints = 10; % Change to whatever you like
    end
    if numberTestPoints == model.N
        testInd = 1:model.N;
    else
        perm = randperm(model.N);
        testInd = perm(1:numberTestPoints);
    end
else
    numberTestPoints = size(Yts{1},1); % 10;
    %perm = randperm(size(Yts{obsMod},1));
    %testInd = perm(1:numberTestPoints);
    if ~exist('testInd')
        testInd = 1:size(Yts{1},1); %%%%  
    end
end

% Here we'll store the predictions: it's a matrix N* x D_z, where N* is the
% number of test points and D_z is the dimensionality of the z points
ZpredMuAll = [];
% Here we'll store the latent points corresponding to the test outputs. If
% the test outputs are totally unobserved, these latent points are found
% via variational optimisation, otherwise (if the test points are actually
% just a subset of the training points), these latent points are found from
% the training latent points already found.
XpredAll = zeros(length(testInd), model.q);  %%%
% The variances of the test latent points
varXpredAll = []; % zeros(length(testInd), model.q);  %%%

if saveAllNN
    % This matrix saves all results, ie is a 1 x N* cell array where each
    % cell contains the top K predictions on z given a y*, where K =
    % numberOfNN.
    ZpredMuAllNN = cell(1, length(testInd));
end
pb = myProgressBar(length(testInd), min(length(testInd),20));

if ~exist('x_star_all', 'var')
    optX = true;
    x_star_all = zeros(length(testInd), model.q);
    varx_star_all = zeros(length(testInd), model.q);
else
    optX = false;
end

for i=1:length(testInd)
    pb = myProgressBar(pb,i);
    curInd = testInd(i);
   % fprintf('# Testing indice number %d ', curInd);
    if testOnTraining
     %   fprintf('taken from the training set\n');
     % y_star is just a copy of a training point y
        y_star = model.comp{obsMod}.y(curInd,:);
     % we don't need to optimise to find the latent point... since y is in
     % the tr. data, then the corresponding x is already found
        x_star = model.comp{obsMod}.vardist.means(curInd,:);
        varx_star = model.comp{obsMod}.vardist.covars(curInd,:);
    else
      %  fprintf('taken from the test set\n');
        y_star = Yts{obsMod}(curInd,:);
        z_star = Yts{infMod}(curInd,:);
        dst = dist2(y_star, model.comp{obsMod}.y);
        [mind, mini(i)] = min(dst);
        
        if ~optX
            x_star = x_star_all(i,:);
            varx_star = varx_star_all(i,:);
        else
            model.comp{obsMod}.mini = mini(i);
            Init(i,:) = model.vardist.means(mini(i),:);
            vardistx = vardistCreate(model.comp{obsMod}.vardist.means(mini(i),:), model.q, 'gaussian');
            vardistx.covars = model.comp{obsMod}.vardist.covars(mini(i),:);
            model.comp{obsMod}.vardistx = vardistx;
            iters = globalOpt.reconstrIters;
            % Find p(X_* | Y_*) which is approximated by q(X_*)
            [x_star, varx_star, modelUpdated] = vargplvmOptimisePoint(model.comp{obsMod}, vardistx, y_star, displayTestOpt, iters);%%%
            x_star_all(i,:) = x_star;
            varx_star_all(i,:) = varx_star;
        end
    end

    % Now we've found a set of datapoints X_* corresponding to the given
    % Y_*. We want to now find the numberOfNN closest nearest neighbours of
    % X_* with the training X but only taking into account the shared
    % dimensions between modalities Y and Z..
   % fprintf('# Finding the %d NN of X_* with the training X based only on the shared dims.\n', numberOfNN);
    [ind, distInd] = nn_class(model.X(:,sharedDims), x_star(:,sharedDims), numberOfNN, 'euclidean');
    
    % Initialise: These matrices hold the numberOfNN predictions for one
    % point z, hence they are numberOfNN x D_z
    ZpredMu = zeros(length(ind), size(model.comp{infMod}.y,2));
    ZpredSigma = zeros(length(ind), size(model.comp{infMod}.y,2));
      
    
    % Find p(y_*|x_*) for every x_* found from the NN
   % fprintf('# Predicting from the NN of X_* ');
    for k=1:numberOfNN
        switch infMethod
            case 0
                % This is the most reliable method (use varx_star in the
                % predictions)
                x_cur = x_star;
            case 1
                x_cur = model.X(ind(k),:);
                x_cur([sharedDims privateDims{obsMod}]) = x_star([sharedDims privateDims{obsMod}]);
            case 2
                x_cur = model.X(ind(k),:);
                x_cur(sharedDims) = x_star(sharedDims); %%% OPTIONAL!!!
            case 3
                %---- Method 3: Combination of predicted x_star and training NN X
                % according to lengthscales (for the shared dims)
                 x_cur = model.X(ind(k),:);
                 s1 = model.comp{obsMod}.kern.comp{1}.inputScales / max(model.comp{obsMod}.kern.comp{1}.inputScales);
                 xcurOrig  = x_cur(sharedDims);
                 s1new = s1/sum(s1);
                 x_cur(sharedDims) = s1new(sharedDims).*x_star(sharedDims) + (1-s1new(sharedDims)).*xcurOrig;
                 x_cur(privateDims{obsMod}) = x_star(privateDims{obsMod});
                %----
            case 4
                x_cur(sharedDims) = x_star(sharedDims);
                x_cur = svargplvmOptimiseDimVar(model.comp{infMod}, x_cur, [privateDims{1} privateDims{2}], 0, 400, 0, 'scg');
            case 5
                % This assumes that there's a GP model which maps 
                % [Xy Xyz] -> Xz, assuming that we do inference for the
                % modality z.
                x_cur = x_star;
                dimsObs = sort([privateDims{obsMod}, sharedDims]);
                dimsInf = privateDims{infMod};
                [x_starGP, x_starGPvar] = gpPosteriorMeanVar(modelGP2, x_cur(:,dimsObs));
                x_cur(:, dimsInf) = x_starGP;
                varx_star(:, dimsInf) = x_starGPvar;
        end
        % Better to use varx_star when x_cur is just taken from x_star
        %ZpredMu(k,:) = vargplvmPosteriorMeanVar(model.comp{infMod}, x_cur);
        ZpredMu(k,:) = vargplvmPosteriorMeanVar(model.comp{infMod}, x_cur, varx_star);%, varx_star(i,:)); % varx_star needed?
    end
    
    % Here we will only hold the 1st (most probable) nearest neighbour
    % found for all test points.
    ZpredMuAll{i} = ZpredMu;%(1,:);
    varXpredAll(i,:) = varx_star;
    XpredAll(i,:) = x_cur;
    
    if saveAllNN
        ZpredMuAllNN{i} = ZpredMu;
    end
        % fprintf('\n\n');
end
fprintf(1, '\n');
