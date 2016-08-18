% VARGPLVMPREDICTUNCERTAIN Extrapolation using an auto-regressive BGPLVM
% DESC
% GENERAL: Given a timeseries x1, x2, ..., xN, we construct a dataset of
% input-output pairs:
% ([x1, ..., xL], x_{L+1}), ([x2, ..., x_{L+1}], x_{L+2}), ... .
% where L is a window size.
% Then, we train a BGPLVM model where the training inputs are fixed and
% have variance 0.
% THIS FUNCTION: Takes a trained model and predicts in the future as
% follows:
% 1. Predict y_{N+1} and get y_{N+1} and vary_{N+1}.
% 2. if augmentModel:
%      2.1: Put y_{N+1} in the input set X as x_{N+1} with variance
%        vary_{N+1}. Pretend this is now a new model trained until point N+1.
%      2.2 if reOptimise:
%        reOptimise the model to get new kernel params.
%    end
% 3. Repeat 2. for Nstar steps (here Nstar is given as a parameter N).
%
%
% Note: In (re)optimisation, one can re-learn X, X_u, both or none.
% Note2: If augmentModel == 1, then even if reOptimiseIters == 0, the model
% will become bigger and bigger so it'll take more and more time to do the
% updates.
%
% ARG model: The trained auto-regressive model
% ARG x_start: The initial input of the future predictions (e.g. if trained
%   until x_n, x_start would be x_{n+1}, or x_test(1,:).
% ARG N: Predict N points ahead
% ARG augmentModel: whether to add predictions in the model and for a new
% "training" model, where "training" refers to what is taken account for
% inference (RECOMMENDED).
% ARG reOptimiseIters: If model is augmented, it can also be reoptimised.
% Give 0 to skip this.
% RETURN Y: The predictive outputs (of size N x d)
% RETURN varY: The predictive variance of the above
% RETURN model: The model (in case it gets updated). The
% model.vardist.covars indices corresponding to newly inserted points
% should be the same as varY.
% RETURN SNR: The SNR of each training procedure (N times) in case
% reOptimiseIters != 0).
%
% NOTE 3: In general, it'd also work if
% a) y (and consequently x) is multivariate and
% b) if instead of predicting e.g. y3 from y1 and y2 we predict (y3 y4) from
% y1 and y2.
% In any case, the above would work with flattening the
% input-output vectors.
% HOWEVER THIS HAS NOT BEEN CAREFULLY TESTED YET (especially b).
%
% NOTE 4: It's not clear that varY should be increasing or decreasing. On
% the one hand, we propagate uncertainty in the inputs. On the other hand,
% the "training set" increases.
%
% SEEALSO: demKstepAhead.m, util_transformSeqToTimeSeries.m,
% util_transformTimeSeriesToSeq.m
function [Y, varY, model, SNR] = vargplvmPredictUncertainK(model, x_start, N, opt)

if nargin < 4,  opt = struct; end
if ~isfield(opt, 'k') || isempty(opt.k)
    k = 1;
else
    k = opt.k;
end
if ~isfield(opt, 'augmentModel') || isempty(opt.augmentModel)
    augmentModel = true; 
else
    augmentModel = opt.augmentModel;
end
if ~isfield(opt, 'augmentInd') || isempty(opt.augmentInd)
    augmentInd = false;
else
    augmentInd = opt.augmentInd;
end
if ~isfield(opt, 'reOptimiseIters') || isempty(opt.reOptimiseIters)
    reOptimiseIters = 20;
else
    reOptimiseIters = opt.reOptimiseIters;
end
if ~isfield(opt, 'trCovarsMult') || isempty(opt.trCovarsMult)
    trCovarsMult = 1;
else
    trCovarsMult = opt.trCovarsMult;
end
if ~isfield(opt, 'varInPred') || isempty(opt.varInPred)
    varInPred = true;
else
    varInPred = opt.varInPred;
end


ws = model.windowSize;
D = model.d;

Y = NaN(N, model.d);
varY = NaN(N, model.d);
curVar = zeros(1, model.q);
curX = x_start;

pb = myProgressBar(N, min(N,50));
if augmentModel && ~isempty(reOptimiseIters) && reOptimiseIters > 0
    SNR = NaN(1,N);
end

m2 = model;
for i=1:N

    %fprintf('%d\n',i); pb.printAll = true; %%% DEBUG
    %%
    if augmentModel && ~mod(i,k)
        m2 = vargplvmUpdateStats(m2, m2.X_u);
        
        if ~isempty(reOptimiseIters) && reOptimiseIters > 0
            %model = vargplvmOptimiseModel(model, 0, 0, {20, reOptimiseIters}, 0);
            m2 = vargplvmOptimise(m2, 0, reOptimiseIters);
            SNR(i) = vargplvmShowSNR(m2,0);
        end
        model = m2;
    end
    if varInPred
        [Ypred, varYpred] = vargplvmPosteriorMeanVar(model, curX, curVar);
    else
        [Ypred, varYpred] = vargplvmPosteriorMeanVar(model, curX);
    end
    
    Y(i,:) = Ypred;
    varY(i,:) = varYpred;
    
    
    
    % ---------- augment model (but not use it unless i == mult of k)
    %j = 0;

    m2.N = m2.N + 1;
    %m2.k = m2.k + 1;
    m2.vardist.means = [m2.vardist.means; curX];
    if trCovarsMult == 1
        m2.vardist.covars = [m2.vardist.covars; curVar];
    else
        m2.vardist.covars = [m2.vardist.covars; trCovarsMult*ones(size(curVar))];
    end
    %m2.vardist.covars = [m2.vardist.covars; curVar+1/m2.beta]; %%%%%%%%%%%%%%%%%%
    
    % Add ind pt if augmentInd is true or if it's not logical but specifies
    % every how many steps we need to add the ind. point.
    if (islogical(augmentInd) && augmentInd)
        m2.X_u = [m2.X_u; curX];
        m2.k = m2.k + size(curX,1);
         pb.printAll = false; %%%%%
    elseif  (~islogical(augmentInd) && ~mod(i,augmentInd))
        X_u = [m2.X_u; curX];
        res = util_checkCloseMatrixElements(X_u, 0.001);
        if ~ismember(size(X_u,1), res)
            % Only add ind. pt if it's not close to another one
            m2.X_u = X_u;
            m2.k = m2.k + size(curX,1);
            %pb.nextSymbol = 'o';
            fprintf('# Adding inuding point...\n');  %%% DEBUG
        else
            fprintf('# Skipping inducing point...\n');  %%% DEBUG
        end
        pb.printAll = true; %%%%%
    end
    m2.X = [m2.X; curX];
    m2.vardist.numData = m2.N;
    m2.vardist.nParams = 2*prod(size(m2.vardist.means));
    m2.vardist.transforms.index = m2.N*m2.q+1:m2.vardist.nParams;
    m2.numParams = m2.numParams + 2*m2.q;
    m2.nParams = m2.numParams;
    
    % normalize y exactly as model.m is normalized
    my = Ypred - repmat(m2.bias,size(Ypred,1),1);
    my = my./repmat(m2.scale,size(Ypred,1),1);
    
    % change the data (by including the new point and taking only the present indices)
    m2.m = [m2.m; my];
    m2.TrYY = sum(sum(m2.m .* m2.m));
    if isfield(m2, 'fixParamIndices')
        m2.fixParamIndices = 1:2*m2.N*m2.q;
    end
    if isfield(m2, 'fixInducing') && m2.fixInducing
        m2.inducingIndices = 1:m2.k;
    end
    %-----------
   % if augmentModel
        pb = myProgressBar(pb,i);
   % end
    
    curX = curX(:,D+1:end);
    curX = [curX Ypred];
    
    curVar = curVar(:, D+1:end);
    curVar = [curVar varYpred];
    
end
fprintf('\n');
if augmentModel && ~isempty(reOptimiseIters) && reOptimiseIters > 0
    fprintf('# SNRs: %s',num2str(SNR));
end





