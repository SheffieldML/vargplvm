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
function [Y, varY, model, SNR] = vargplvmPredictUncertain(model, x_start, N, augmentModel, reOptimiseIters)

if nargin < 4 || isempty(augmentModel), augmentModel = true; end
if nargin < 5 || isempty(reOptimiseIters), reOptimiseIters = 20; end

ws = model.windowSize;
D = model.d;

Y = NaN(N, model.d);
varY = NaN(N, model.d);
curVar = zeros(1, model.q);
curX = x_start;

if augmentModel %&& ~isempty(reOptimiseIters) && reOptimiseIters > 0
    pb = myProgressBar(N, min(N,50));
    SNR = NaN(1,N);
end

for i=1:N    
    [Ypred, varYpred] = vargplvmPosteriorMeanVar(model, curX, curVar);
    Y(i,:) = Ypred;
    varY(i,:) = varYpred;
    
    %%
    if augmentModel
        %j = 0;
        model.N = model.N + 1;
        model.k = model.k + 1;
        model.vardist.means = [model.vardist.means; curX];
        model.vardist.covars = [model.vardist.covars; curVar];
        model.X_u = [model.X_u; curX];
        model.X = [model.X; curX];
        model.vardist.numData = model.N;
        model.vardist.nParams = 2*prod(size(model.vardist.means));
        model.vardist.transforms.index = model.N*model.q+1:model.vardist.nParams;
        model.numParams = model.numParams + 2*model.q;
        model.nParams = model.numParams;
        
        % normalize y exactly as model.m is normalized
        my = Ypred - repmat(model.bias,size(Ypred,1),1);
        my = my./repmat(model.scale,size(Ypred,1),1);
        
        % change the data (by including the new point and taking only the present indices)
        model.m = [model.m; my];
        model.TrYY = sum(sum(model.m .* model.m));
        if isfield(model, 'fixParamIndices')
            model.fixParamIndices = 1:2*model.N*model.q;
        end
        if isfield(model, 'fixInducing') && model.fixInducing
            model.inducingIndices = 1:model.k;
        end
        
        model = vargplvmUpdateStats(model, model.X_u);
        
        if ~isempty(reOptimiseIters) && reOptimiseIters > 0
            %model = vargplvmOptimiseModel(model, 0, 0, {20, reOptimiseIters}, 0);
            model = vargplvmOptimise(model, 0, reOptimiseIters);
            SNR(i) = vargplvmShowSNR(model,0);
        end
        pb = myProgressBar(pb,i);
    end
    
    curX = curX(:,model.d+1:end);
    curX = [curX Ypred];
    
    curVar = curVar(:, model.d+1:end);
    curVar = [curVar varYpred];
    
end
if augmentModel && ~isempty(reOptimiseIters) && reOptimiseIters > 0
    fprintf('# SNRs: %s',num2str(SNR));
end





