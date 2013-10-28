function [X, varX, model, grChek] = svargplvmOptimisePoint(model, vardistx, y, display, iters, varargin)

% SVARGPLVMOPTIMISEPOINT Optimise the postion of one or more latent points
% given observations from 1 or more modalities
% FORMAT
% DESC optimises the location of a group of points in latent space
% given an initialisation and the corresponding observed data point.
% ARG model : the MRD model for which the point will be optimised.
% ARG vardistx : the initialisation of the points in the latent space.
% ARG y : the observed data points for which the latent points are to
% be optimised.
% ARG display : whether or not to display the iterations of the
% optimisation (default: true)
% ARG iters : maximum number of iterations for the optimisation
% (default 2000).
% RETURN x : the optimised means in the latent space.
% RETURN varx : the optimised variances in the latent space.
% RETURN model: the model which is augmented to also include the new test
% points and the quantities that change because of these, as there is
% coupling in the dynamics case.

%
% COPYRIGHT :  Andreas Damianou, 2013
%
% SEEALSO : vargplvmOptimisePoint

% VARGPLVM

grChek = [];

if nargin < 5 || isempty(iters), iters = 2000; end
if nargin < 4 || isempty(display), display = true; end

options = optOptions;
if display
    options(1) = 1;
    % options(9) = 1; % gradchek
end
options(14) = iters;


if isfield(model, 'optimiser')
    optim = str2func(model.optimiser);
else
    optim = str2func('scg');
end


if isfield(model, 'dynamics') && ~isempty(model.dynamics)
    % TODO
    error('Not implemented for dynamical case')
else
    x = vardistExtractParam(vardistx);
end

if ~(isfield(model, 'dynamics') && ~isempty(model.dynamics))
    model.DgtN_test = true;
    for i = model.testModalities
        % Perform precomputations which will result in faster execution. This
        % happens in two stages:
        % a) Precomputing constants only once here, instead in every call of
        % the objective and grad
        % b) Replace model.m with a reduced rank representation, since it only
        % appears in the expression model.m * model.m'. This is the same tricks
        % used in vargplvmCreate for DgtN flag.
        % TODO: Do the same for dynamics case
        if isfield(model.comp{i}, 'DgtN') && model.comp{i}.DgtN
            mOrig{i} = model.comp{i}.m;
        end
        model.comp{i}.testPrecomp.indexMissing = find(isnan(y{i}(1,:)));
        indexPresent = setdiff(1:model.comp{i}.d, model.comp{i}.testPrecomp.indexMissing );
        y{i} = y{i}(:,indexPresent);
        P = model.comp{i}.P1 * (model.comp{i}.Psi1' * model.comp{i}.m(:,model.comp{i}.testPrecomp.indexMissing));
        model.comp{i}.testPrecomp.TrPP = sum(sum(P .* P));
        model.comp{i}.testPrecomp.TrYY = sum(sum(model.comp{i}.m(:,model.comp{i}.testPrecomp.indexMissing) .* model.comp{i}.m(:,model.comp{i}.testPrecomp.indexMissing)));
        y{i} = y{i} - repmat(model.comp{i}.bias(indexPresent),size(y{i},1),1);
        y{i} = y{i}./repmat(model.comp{i}.scale(indexPresent),size(y{i},1),1);
        mPres = model.comp{i}.m(:, indexPresent);
        mPres = [mPres; y{i}];
        YYT = mPres * mPres'; % NxN
        [U S V]=svd(YYT);
        if isfield(model.comp{i}, 'DgtN') && model.comp{i}.DgtN
            model.comp{i}.testPrecomp.mReduced=U*sqrt(abs(S));
        end
        model.comp{i}.testPrecomp.TrYY2 = sum(diag(YYT)); % scalar
        
        % For grads
        mPresGrad = [y{i}; model.comp{i}.m(:, indexPresent)];
        YYT = mPresGrad * mPresGrad'; % NxN
        [U S V]=svd(YYT);
        if isfield(model.comp{i}, 'DgtN') && model.comp{i}.DgtN
            model.comp{i}.testPrecomp.mReducedGrad=U*sqrt(abs(S));
        end
        model.comp{i}.testPrecomp.mY = mPresGrad*y{i}';
        if isfield(model.comp{i}, 'DgtN') && model.comp{i}.DgtN
            model.comp{i}.m = []; % Less overhead in passing arguments (pass by value)
        end
    end
else
    model.comp.DgtN_test = false;
end



if length(varargin) == 2
    if strcmp(varargin{1}, 'gradcheck')
        assert(islogical(varargin{2}));
        %options(9) = varargin{2};
        doGradchek = varargin{2};
        if doGradchek
            [gradient, delta] = feval('gradchek', vardistExtractParam(vardistx), @svargplvmPointObjective, @svargplvmPointGradient, model, y);
            deltaf = gradient - delta;
            d=norm(deltaf - gradient)/norm(gradient + deltaf); %%
            d1=norm(deltaf - gradient,1)/norm(gradient + deltaf,1); %%
            grRatio = sum(abs(gradient ./ deltaf)) / length(deltaf);
            fprintf(1,' Norm1 difference: %d\n Norm2 difference: %d\n Ratio: %d\n',d1,d, grRatio);
            grChek = {delta, gradient, deltaf, d, d1};
        else
            grChek = [];
        end
    end
end


if iters > 0
    if strcmp(func2str(optim), 'optimiMinimize')
        % Carl Rasmussen's minimize function
        x = optim('vargplvmPointObjectiveGradient', x, options, model, y);
    else
        % NETLAB style optimization.
        x = optim('svargplvmPointObjective', x,  options, ...
            'svargplvmPointGradient', model, y);
    end
end

if isfield(model, 'dynamics') && ~isempty(model.dynamics)
    if isfield(model.dynamics, 'reoptimise') && model.dynamics.reoptimise   %%% RE-OPT-CODE-NEW
        % Expand and update the model.                                      %%% RE-OPT-CODE-NEW
        % In this case, x is [mu_bar lambda theta_t X_u]                    %%% RE-OPT-CODE-NEW
        [vardistx, model] = vargplvmPartExpand(model, x, 1);                %%% RE-OPT-CODE-NEW
    else %%% RE-OPT-CODE-NEW
        % now separate the variational disribution into the training part and the
        % testing part and update the original training model (only with the new training
        % variational distribution) and the test variational distribution
        % this is doing the expand
        x = reshape(x, vardist.numData, model.dynamics.q*2);
        xtrain = x(1:model.N,:);
        xtest = x(model.N+1:end,:);
        model.dynamics.vardist = vardistExpandParam(model.dynamics.vardist, xtrain);
        vardistx = vardistExpandParam(model.vardistx, xtest);
    end %%% RE-OPT-CODE-NEW
else
    vardistx = vardistExpandParam(vardistx,x);
end

X = vardistx.means;
varX = vardistx.covars;

if ~(isfield(model, 'dynamics') && ~isempty(model.dynamics))
    model = rmfield(model, 'DgtN_test');
    for i=model.testModalities
        model.comp{i} = rmfield(model.comp{i}, 'testPrecomp');
        if isfield(model.comp{i}, 'DgtN') && model.comp{i}.DgtN
            model.comp{i}.m = mOrig{i};
        end
    end
end
