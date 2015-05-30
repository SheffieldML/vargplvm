function [x_star_all, varx_star_all, mini] = ...
    vargplvmPredictLatent(model, Yts, testInd, batchMode, iters, displayTestOpt, mini, initQx, parallel)

% VARGPLVMPREDICTLATENT Given output observation(s) compute the approximate
% posterior
% DESC Given an output observed test point y* (or set of test points Y*)
% compute the approximate posterior p(X* | Y*) for the corresponding latent
% points.
%
% ARG model: The trained vargplvm model
% ARG Yts: The set of outputs observations
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
% computed here in the default way.
% ARG initQx: Direct initialisation (instead of mini). Has to contain
% .means and .covars
%
% SEEALSO:
%
% COPYRIGHT : Andreas C. Damianou, 2013
%
% VARGPLVM

if nargin < 2 || isempty(model) || isempty(Yts), error('At least 2 arguments needed!'); end
if nargin < 3 || isempty(testInd), testInd = 1:size(Yts,1); end
if nargin < 4 || isempty(batchMode), batchMode = false; end
if nargin < 5 || isempty(iters), iters = model.globalOpt.reconstrIters; end
if nargin < 6 || isempty(displayTestOpt), displayTestOpt = false; end
if nargin < 7, mini = []; end
if nargin < 8, initQx = []; end
if nargin < 9, parallel = []; end

try
    pool_open = matlabpool('size')>0;
catch e
    pool_open = 0;
end
if pool_open && ~isempty(parallel) && size(Yts,1) > 3 && ~batchMode
    parallel = parallel;
    warning('Parallel predictions!')
else
    parallel = false;
end

miniGiven = ~isempty(mini);
initQxGiven = ~isempty(initQx);

if initQxGiven && miniGiven
    error('Both initial training indices AND initial X was given!')
end

if batchMode
    y_star = Yts(testInd,:);
    if ~initQxGiven
        if ~miniGiven
            mini = NaN(1,length(testInd));
            for i=1:length(testInd)
                dst = dist2(y_star(i,:), model.y);
                [~, mini(i)] = min(dst);
            end
        end
        model.mini = mini;
        vardistx = vardistCreate(model.vardist.means(mini, :), model.q, 'gaussian');
        vardistx.covars = model.vardist.covars(mini, :);
    else
        vardistx = vardistCreate(initQx.means, model.q, 'gaussian');
        vardistx.covars = initQx.covars;
    end
    model.vardistx = vardistx;
    [x_star_all, varx_star_all] = vargplvmOptimisePoint(model, vardistx, y_star, true, iters);
else
    if parallel
        x_star_allCell = cell(length(testInd), 1);
        varx_star_allCell = cell(length(testInd), 1);
        if ~initQxGiven &&  ~miniGiven
            mini = NaN(1,length(testInd));
        end
        for i=1:length(testInd)
             curInd = testInd(i);
             YYts{i} = Yts(curInd,:);
             if initQxGiven
                initQxMeans{i} = initQx.means(i,:);
                initQxCovars{i} = initQx.covars(i,:);
             else
                 initQxMeans{i} = [];
                 initQxCovars{i} = [];
             end
             miniCell{i} = mini(i);
        end
         Ytr = model.y;
        pb = myProgressBar(length(testInd), min(length(testInd),20));
        for i=1:length(testInd)
            pbb{i} = pb;
        end
        parfor i=1:length(testInd)
            pbb{i} = myProgressBar(pbb{i},i);
            y_star =YYts{i};
            if ~initQxGiven
                if ~miniGiven
                    dst = dist2(y_star, Ytr);
                    [~, miniCur] = min(dst);
                    miniCell{i} = miniCur;
                else
                    miniCur = miniCell{i};
                end
                modelCur = model;
                modelCur.mini = miniCur;
                vardistx = vardistCreate(modelCur.vardist.means(miniCur,:), modelCur.q, 'gaussian');
                vardistx.covars = modelCur.vardist.covars(miniCur,:);
            else
                vardistx = vardistCreate(initQxMeans{i}, modelCur.q, 'gaussian');
                vardistx.covars = initQxCovars{i};
            end
            modelCur.vardistx = vardistx;
            
            % Find p(X_* | Y_*) which is approximated by q(X_*)
            [x_star, varx_star] = vargplvmOptimisePoint(modelCur, vardistx, y_star, displayTestOpt, iters);
            x_star_allCell{i} = x_star;
            varx_star_allCell{i} = varx_star;
        end
        for i=1:length(testInd)
            x_star_all(i,:) = x_star_allCell{i};
            varx_star_all(i,:) = varx_star_allCell{i};
            mini(i) = miniCell{i};
        end
        fprintf(1, '\n');
%         %.. TODO
%         x_star_all = zeros(length(testInd), model.q);
%         varx_star_all = zeros(length(testInd), model.q);
%         if ~initQxGiven &&  ~miniGiven
%             mini = NaN(1,length(testInd));
%         end
%         pb = myProgressBar(length(testInd), min(length(testInd),20));
%         YYts = cell(length(testInd),1);
%         initQxMeans = cell(length(testInd),1);
%         initQxCovars = cell(length(testInd),1);
%         for i=1:length(testInd)
%             curInd = testInd(i);
%             YYts{i} = Yts(curInd,:);
%             initQxMeans{i} = initQx.means(i,:);
%             initQxCovars{i} = initQx.covars(i,:);
%         end
%         parfor i=1:length(testInd)
%             pb = myProgressBar(pb,i);
%             %curInd = testInd(i);
%             y_star = YYts{i};
%             modelTmp = model;
%             if ~initQxGiven
%                 if ~miniGiven
%                     dst = dist2(y_star, model.y);
%                     [~, mini(i)] = min(dst);
%                 end
%                 modelTmp.mini = mini(i);
%                 vardistx = vardistCreate(modelTmp.vardist.means(mini(i),:), modelTmp.q, 'gaussian');
%                 vardistx.covars = modelTmp.vardist.covars(mini(i),:);
%             else
%                 vardistx = vardistCreate(initQxMeans{i}, modelTmp.q, 'gaussian');
%                 vardistx.covars = initQxCovars{i};
%             end
%             modelTmp.vardistx = vardistx;
%             
%             % Find p(X_* | Y_*) which is approximated by q(X_*)
%             [x_star, varx_star] = vargplvmOptimisePoint(modelTmp, vardistx, y_star, displayTestOpt, iters);
%             x_star_all(i,:) = x_star;
%             varx_star_all(i,:) = varx_star;
%         end
%         fprintf(1, '\n');     
    else
        x_star_all = zeros(length(testInd), model.q);
        varx_star_all = zeros(length(testInd), model.q);
        if ~initQxGiven &&  ~miniGiven
            mini = NaN(1,length(testInd));
        end
        pb = myProgressBar(length(testInd), min(length(testInd),20));
        for i=1:length(testInd)
            pb = myProgressBar(pb,i);
            curInd = testInd(i);
            y_star = Yts(curInd,:);
            if ~initQxGiven
                if ~miniGiven
                    dst = dist2(y_star, model.y);
                    [~, mini(i)] = min(dst);
                end
                model.mini = mini(i);
                vardistx = vardistCreate(model.vardist.means(mini(i),:), model.q, 'gaussian');
                vardistx.covars = model.vardist.covars(mini(i),:);
            else
                vardistx = vardistCreate(initQx.means(i,:), model.q, 'gaussian');
                vardistx.covars = initQx.covars(i,:);
            end
            model.vardistx = vardistx;
            
            % Find p(X_* | Y_*) which is approximated by q(X_*)
            [x_star, varx_star] = vargplvmOptimisePoint(model, vardistx, y_star, displayTestOpt, iters);
            x_star_all(i,:) = x_star;
            varx_star_all(i,:) = varx_star;
        end
        fprintf(1, '\n');
    end
end