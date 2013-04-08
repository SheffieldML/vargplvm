function model = svargplvmOptimise(model, display, iters, varargin)

% SHEFFIELDMLOPTIMISE Optimise the SVARGPLVM.
% FORMAT
% DESC takes a given shared variational GP-LVM model structure and optimises with
% respect to parameters and latent positions. 
% ARG model : the model to be optimised.
% ARG display : flag dictating whether or not to display
% optimisation progress (set to greater than zero) (default value 1). 
% ARG iters : number of iterations to run the optimiser
% for (default value 2000).
% RETURN model : the optimised model.
%
% SEEALSO : svargplvmModelCreate, svargplvmLogLikelihood,
% svargplvmLogLikeGradients, svargplvmObjective, svargplvmGradient,
% svargplvmOptimiseModel
% 
% COPYRIGHT: Andreas C. Damianou, 2011

% SHEFFIELDML


if nargin < 3
  iters = 2000;
  if nargin < 2
    display = 1;
  end
end

options = optOptions;
if length(varargin) == 2
    if strcmp(varargin{1}, 'gradcheck')
        assert(islogical(varargin{2}));
        options(9) = varargin{2};
        if options(9)
            [params, names] = svargplvmExtractParam(model);
            for i=1:length(names)
                fprintf('%d\t%s\n', i, names{i});
            end
        else
            params = svargplvmExtractParam(model);
        end
    end
else
    params = svargplvmExtractParam(model);
end

%params = svargplvmExtractParam(model);

options(2) = 0.1*options(2); 
options(3) = 0.1*options(3);

if display
  options(1) = 1;
  if length(params) <= 100
    options(9) = 1; % gradchek
  end
end
options(14) = iters;

if isfield(model, 'optimiser')
  optim = str2func(model.optimiser);
else
  optim = str2func('scg');
end


if strcmp(func2str(optim), 'optimiMinimize')
  % Carl Rasmussen's minimize function 
  params = optim('svargplvmObjectiveGradient', params, options, model);
elseif strcmp(func2str(optim), 'scg2')
   % NETLAB style optimization with a slight modification so that an
   % objectiveGradient can be used where applicable, in order to re-use
   % precomputed quantities.
   params = optim('svargplvmObjectiveGradient', params,  options, 'svargplvmGradient', model);
else
   % NETLAB style optimization.
   params = optim('svargplvmObjective', params,  options,  'svargplvmGradient', model);
end

model = svargplvmExpandParam(model, params);

