function [x_star_all, varx_star_all, mini] = ...
    svargplvmPredictLatent(model, Yts, testInd, obsMod, batchMode, iters, displayTestOpt, mini)

% SVARGPLVMPREDICTLATENT
%
% DESC Given one ore more output observations find the corresponding latent
% point from the approximate posterior.
% Similar to vargplvmPredictLatent but this function accepts more than
% one given modalities (if only one given it should be identical to
% vargplvmPredictLatent).
%
% ARG model: The trained svargplvm (MRD) model
% ARG Yts: The set of outputs observations. It's a cell array with one cell
% being the matrix of observations for one modality
% ARG testInd: The indices of Yts to do inference for (default: all)
% ARG batchMode: If set to true, then if N* points are given in Yts, only
% one optimisation will be done in the N*-dimensional space. Otherwise, N*
% 1-dimensional optimisations will be performed. Batch mode is faster but
% possibly less accurate.
% ARG iters: the number of test iterations for the optimisation(s)
% ARG displayTestOpt: If set to true, then the test iterations will be
% displayed
% ARG mini: The indices of the points in the test set for which the test
% latent points are initialised with, ie x*(i,:) will be initialised with
% Xtrain(mini(i), :). If this vector is not given then this vector is
% computed here in the default way. Alternatively, if mini is a cell array
% of length 2, then mini{1} represents the initial x* and mini{2}
% represents their initial associated variances.
%
% SEEALSO : vargplvmPredictLatent
%
% COPYRIGHT : Andreas C. Damianou, 2013
%
% VARGPLVM

if nargin < 2 || isempty(model) || isempty(Yts), error('At least 2 arguments needed!'); end
if nargin < 3 || isempty(testInd), testInd = 1:size(Yts{1},1); end
if nargin < 4 || isempty(obsMod), obsMod = 1:length(Yts); end
if nargin < 5 || isempty(batchMode), batchMode = false; end
if nargin < 6 || isempty(iters), iters = model.globalOpt.reconstrIters; end
if nargin < 7 || isempty(displayTestOpt), displayTestOpt = false; end
if nargin < 8, mini = []; end

miniGiven = ~isempty(mini);

for jj = obsMod
    Yts{jj} = Yts{jj}(testInd, :);
end

if ~isfield(model, 'testModalities')
    model.testModalities = obsMod;
end

if batchMode
    if ~miniGiven
        y_star = [];
        tmpYtr = [];
        for jj = obsMod
            y_star = [y_star Yts{jj}];
            tmpYtr = [tmpYtr model.comp{jj}.y];
            
        end
        mini = NaN(1,length(testInd));
        for i=1:length(mini)
            dst = dist2(y_star(i,:), tmpYtr);
            [~, mini(i)] = min(dst);
        end
    end
    if iscell(mini)
        vardistx = vardistCreate(mini{1}, model.q, 'gaussian');
        vardistx.covars = mini{2};
    else
        model.mini = mini;
        vardistx = vardistCreate(model.vardist.means(mini, :), model.q, 'gaussian');
        vardistx.covars = model.vardist.covars(mini, :);
    end
    model.vardistx = vardistx;
    [x_star_all, varx_star_all] = svargplvmOptimisePoint(model, vardistx, Yts, true, iters);
else
    x_star_all = zeros(length(testInd), model.q);
    varx_star_all = zeros(length(testInd), model.q);
    if ~miniGiven
        mini = NaN(1,length(testInd));
    end
    pb = myProgressBar(length(testInd), min(length(testInd),20));
    for i=1:length(testInd)
        pb = myProgressBar(pb,i);  
        if ~miniGiven
            tmpYtr = [];
            y_star = [];
            for jj=obsMod
                y_star = [y_star Yts{jj}(i, :)];
                tmpYtr = [tmpYtr model.comp{jj}.y];
            end
            dst = dist2(y_star, tmpYtr);
            [~, mini(i)] = min(dst);
        end
        if iscell(mini)
            vardistx = vardistCreate(mini{1}(i,:), model.q, 'gaussian');
            vardistx.covars = mini{2}(i,:);
        else
            model.mini = mini(i);
            vardistx = vardistCreate(model.vardist.means(mini(i),:), model.q, 'gaussian');
            vardistx.covars = model.vardist.covars(mini(i),:);
        end
        model.vardistx = vardistx;
        yts_cur = cell(1, length(Yts));
        for jj = obsMod
            yts_cur{jj} = Yts{jj}(i,:);
        end
        % Find p(X_* | Y_*) which is approximated by q(X_*)
        [x_star, varx_star] = svargplvmOptimisePoint(model, vardistx, yts_cur, displayTestOpt, iters);
        x_star_all(i,:) = x_star;
        varx_star_all(i,:) = varx_star;
    end
    fprintf(1, '\n');
end
