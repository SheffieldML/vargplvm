function [ZpredMuAll, testInd, XpredAll, varXpredAll, indNN, ZpredSigmaAll] = ...
    svargplvmPredictionsFunc( ...
        model, testOnTraining, x_star_all, varx_star_all, ...
        obsMod, infMod, testInd, numberOfNN, infMethod, privSharedDims, nnToKeep, modelGP, returnVar)

% SVARGPLVMPREDICTIONSFUNC Identical to svargplvmPredictions but in
% function form
%
% FORMAT: 
%function [ZpredMuAll, testInd, XpredAll, varXpredAll] = ...
%    svargplvmPredictionsFunc( ...
%        model, testOnTraining, x_star_all, varx_star_all,
%        obsMod, infMod, testInd, numberOfNN, infMethod, privSharedDims, nnToKeep)
% 
% ARG model: The svargplvm model
% ARG testOnTraining: true or false, whether to attempt reconstruction on
% training data or predict test data
% ARG x_star_all:
% ARG varx_star_all: 
% ARG obsMod: Yts{infMod} given Yts{obsMod} (Yts is the test set, not given
% as input to this function)
% ARG infMod: Yts{infMod} given Yts{obsMod} 
% ARG testInd: Alternatively to the above, give the indices to
% reconstruct/predict
% ARG numberOfNN: How many points in the other modality we want to recover /
% generate per single test output y*
% ARG infMethod: a number 0-5, different inference methods to find z*|x*
% given x*|y*
% ARG privSharedDims: optional private/shared dimensions of modalities
% ARG nnToKeep: Which NN indice to keep (useful when testing on training
% data). Can be a single index (scalar) or set of indexes (vector)
% ARG returnVar:Whether to return the predictive variance (slower) or not
%
% SEEALSO : demSvargplvmGeneric, demClassification, demClassification3,
% demClassificationGeneral
%
% COPYRIGHT : Andreas C. Damianou, 2013
%
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

if nargin < 1 || isempty(model), error('At least one argument needed!'); end
if nargin < 2 || isempty(testOnTraining), testOnTraining = false; end
if ~testOnTraining && (nargin < 4 || isempty(x_star_all) || isempty(varx_star_all))
    error('Missing latent predictions! Use svargplvmPredictLatent first!');
end
if nargin < 5 || isempty(obsMod), obsMod = 1; end
if nargin < 6 || isempty(infMod), infMod = 2; end
if nargin < 7 || isempty(testInd)
    if testOnTraining
        testInd = 1:model.N;
    else
        testInd = 1:size(x_star_all,1);
    end
end
if nargin < 8 || isempty(numberOfNN)
    if nargin < 11 || isempty(nnToKeep)
        numberOfNN = 1; 
    else
        numberOfNN = max(nnToKeep);
    end
end
if nargin < 9 || isempty(infMethod), infMethod = 1; end
if nargin < 10, privSharedDims = {}; end
if nargin < 11 || isempty(nnToKeep), nnToKeep = 1; end
if nargin < 12, modelGP = []; end
if nargin < 13 ||isempty(returnVar), returnVar = false; end

% TODO:
% Here replace model.comp{i}.m with model.comp{i}.mOrig if DgtN is
% active...

for i=1:model.numModels
    if strcmp(model.comp{i}.kern.type, 'rbfardjit')
        if isfield(model.comp{i}.kern, 'comp')
            assert(sum(abs(model.comp{i}.kern.inputScales - model.comp{i}.kern.comp{1}.inputScales)) == 0);
        else
            model.comp{i}.kern.comp{1}.inputScales = model.comp{i}.kern.inputScales;
        end
    end
end

if isempty(privSharedDims)
    [sharedDims, privateDims] = svargplvmFindSharedDims(model,[],[],{obsMod infMod});
else
    sharedDims = privSharedDims{1};
    privateDims = privSharedDims{2};
end


% Here we'll store the latent points corresponding to the test outputs. If
% the test outputs are totally unobserved, these latent points are found
% via variational optimisation, otherwise (if the test points are actually
% just a subset of the training points), these latent points are found from
% the training latent points already found.
if nargout > 2
    XpredAll = zeros(length(testInd), model.q);
end
if nargout > 3
    % The variances of the test latent points
    varXpredAll = zeros(length(testInd), model.q);
end

% Here we'll store the predictions: it's a matrix N* x D_z, where N* is the
% number of test points and D_z is the dimensionality of the z points
ZpredMuAll = cell(1, length(testInd));

if nargout > 4
    indNN = NaN(length(testInd), length(nnToKeep));
end

pb = myProgressBar(length(testInd), min(length(testInd),20));
for i=1:length(testInd)
    pb = myProgressBar(pb,i);
    curInd = testInd(i);
    if testOnTraining
        x_star = model.vardist.means(curInd,:);
        varx_star = model.vardist.covars(curInd,:);
    else
        x_star = x_star_all(i,:);
        varx_star = varx_star_all(i,:);
    end
    
    % Now we've found a set of datapoints X_* corresponding to the given
    % Y_*. We want to now find the numberOfNN closest nearest neighbours of
    % X_* with the training X but only taking into account the shared
    % dimensions between modalities Y and Z..
   % fprintf('# Finding the %d NN of X_* with the training X based only on the shared dims.\n', numberOfNN);
    [ind, ~] = nn_class(model.X(:,sharedDims), x_star(:,sharedDims), numberOfNN, 'euclidean');
    
    if nargout > 4
        indNN(i, :) = ind(nnToKeep);
    end
    
    % Initialise: These matrices hold the numberOfNN predictions for one
    % point z, hence they are numb    erOfNN x D_z
    ZpredMu = zeros(length(ind), size(model.comp{infMod}.y,2));
    if returnVar
        ZpredSigma = zeros(length(ind), size(model.comp{infMod}.y,2));
    end
      
    allPrivObs = cell2mat(privateDims(obsMod(:)));
                    
    % Find p(y_*|x_*) for every x_* found from the NN
   % fprintf('# Predicting from the NN of X_* ');
    for k=1:numberOfNN
        switch infMethod
            case 0
                % This is the most reliable method (use varx_star in the
                % predictions)
                x_cur = x_star;
            case 1
                % Replace the private dims of the inf. model with the NN
                x_cur = model.X(ind(k),:);
                %x_cur([sharedDims privateDims{obsMod}]) = x_star([sharedDims privateDims{obsMod}]);
                x_cur([sharedDims allPrivObs]) = x_star([sharedDims allPrivObs]);
            case 2
                % Replace the private dims of the inf. model with the NN
                % ones AS WELL AS the private of the observed model
                x_cur = model.X(ind(k),:);
                x_cur(sharedDims) = x_star(sharedDims); %%% OPTIONAL!!!
            case 3
                %---- Method 3: Combination of predicted x_star and training NN X
                % according to lengthscales (for the shared dims)
                 x_cur = model.X(ind(k),:);
                 s1 = model.comp{obsMod(1)}.kern.comp{1}.inputScales / max(model.comp{obsMod(1)}.kern.comp{1}.inputScales);
                 for jj=2:length(obsMod)
                    s1 = (s1 + model.comp{obsMod(jj)}.kern.comp{1}.inputScales / max(model.comp{obsMod(jj)}.kern.comp{1}.inputScales)) / length(obsMod);
                 end
                 xcurOrig  = x_cur(sharedDims);
                 s1new = s1/sum(s1);
                 x_cur(sharedDims) = s1new(sharedDims).*x_star(sharedDims) + (1-s1new(sharedDims)).*xcurOrig;
                 x_cur(allPrivObs) = x_star(allPrivObs);
                %----
            case 4
                x_cur(sharedDims) = x_star(sharedDims);
                x_cur = svargplvmOptimiseDimVar(model.comp{infMod}, x_cur, [privateDims{1} privateDims{2}], 0, 400, 0, 'scg');
            case 5
                % This assumes that there's a GP model which maps 
                % [Xy Xyz] -> Xz, assuming that we do inference for the
                % modality z
                x_cur = x_star;
                dimsObs = sort(unique([allPrivObs, sharedDims])); % TODO if obsMod has > 1 modalities
                dimsInf = privateDims{infMod};
                [x_starGP, x_starGPvar] = gpPosteriorMeanVar(modelGP, x_cur(:,dimsObs));
                x_cur(:, dimsInf) = x_starGP;
                varx_star(:, dimsInf) = x_starGPvar;
            case 6
                % TODO
                % Xz -> [Xy Xyz] % GPLVM
        end
        % Better to use varx_star when x_cur is just taken from x_star
        %ZpredMu(k,:) = vargplvmPosteriorMeanVar(model.comp{infMod}, x_cur);
        if returnVar
            [ZpredMu(k,:) ZpredSigma(k,:)] = vargplvmPosteriorMeanVar(model.comp{infMod}, x_cur, varx_star);
        else
            ZpredMu(k,:) = vargplvmPosteriorMeanVar(model.comp{infMod}, x_cur, varx_star);%, varx_star(i,:)); 
        end
    end
    
    % Here we will only hold the nearest neighbours found for all test points.
    ZpredMuAll{i} = ZpredMu(nnToKeep, :);
    if returnVar
        ZpredSigmaAll{i} = ZpredSigma(nnToKeep,:);
    end
    
    if nargout > 2
        XpredAll(i,:) = x_cur;
    end
    if nargout > 3
        varXpredAll(i,:) = varx_star;
    end
end
fprintf('\n')