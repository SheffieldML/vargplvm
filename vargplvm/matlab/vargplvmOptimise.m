function [model, grChek] = vargplvmOptimise(model, display, iters, varargin)

% VARGPLVMOPTIMISE Optimise the VARGPLVM.
% FORMAT
% DESC takes a given GP-LVM model structure and optimises with
% respect to parameters and latent positions.
% ARG model : the model to be optimised.
% ARG display : flag dictating whether or not to display
% optimisation progress (set to greater than zero) (default value 1).
% ARG iters : number of iterations to run the optimiser
% for (default value 2000).
% RETURN model : the optimised model.
%
% SEEALSO : vargplvmCreate, vargplvmLogLikelihood,
% vargplvmLogLikeGradients, vargplvmObjective, vargplvmGradient
%
% COPYRIGHT : Michalis K. Titsias, 2009
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% VARGPLVM


if nargin < 3
    iters = 2000;
    if nargin < 2
        display = 1;
    end
end

options = optOptions;
params = modelExtractParam(model);
if length(varargin) == 2
    if strcmp(varargin{1}, 'gradcheck')
        assert(islogical(varargin{2}));
        %options(9) = varargin{2};
        doGradchek = varargin{2};
        if doGradchek
            [gradient, delta] = feval('gradchek', params, @vargplvmObjective, @vargplvmGradient, model);
            deltaf = gradient - delta;
            d=norm(deltaf - gradient)/norm(gradient + deltaf);
            d1=norm(deltaf - gradient,1)/norm(gradient + deltaf,1);
            grRatio = sum(abs(gradient ./ deltaf)) / length(deltaf);
            fprintf(1,' Norm1 difference: %d\n Norm2 difference: %d\n Ratio: %d\n',d1,d, grRatio);
            grChek = {delta, gradient, deltaf, d, d1};
        else
            grChek = [];
        end
    end
end

options(2) = 0.1*options(2);
options(3) = 0.1*options(3);

if display
    options(1) = 1;
 %   if length(params) <= 100
 %       options(9) = 1;
 %   end
end
options(14) = iters;

if iters > 0
    if isfield(model, 'optimiser') && ~isa(model.optimiser, 'function_handle')
        if isfield(model, 'optimiser')
            optim = str2func(model.optimiser);
        else
            optim = str2func('scg');
        end
        
        if strcmp(func2str(optim), 'optimiMinimize')
            % Carl Rasmussen's minimize function
            params = optim('vargplvmObjectiveGradient', params, options, model);
        elseif strcmp(func2str(optim), 'scg2')
            % NETLAB style optimization with a slight modification so that an
            % objectiveGradient can be used where applicable, in order to re-use
            % precomputed quantities.
            params = optim('vargplvmObjectiveGradient', params,  options,  'vargplvmGradient', model);
        else
            % NETLAB style optimization.
            params = optim('vargplvmObjective', params,  options,  'vargplvmGradient', model);
        end
    elseif isfield(model, 'optimiser') && isa(model.optimiser, 'function_handle')
        f = fcnchk(model.optimiser);
        params = f(model);
    else
        error('vargplvmOptimise: Invalid optimiser setting.');
    end
    %model = vargplvmExpandParam(model, params);
    model = modelExpandParam(model, params);
    
    % Check if SNR of the optimised model is reasonable (ortherwise a
    % bad local minimum might have been found)
    if isfield(model, 'throwSNRError')
        svargplvmCheckSNR({vargplvmShowSNR(model)}, [], [], model.throwSNRError);
    else
        svargplvmCheckSNR({vargplvmShowSNR(model)});
    end
end