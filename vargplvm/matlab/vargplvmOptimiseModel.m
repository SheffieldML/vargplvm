function model = vargplvmOptimiseModel(model, varargin)

% VARGPLVMOPTIMISEMODEL Optimise a varglpvm model (wrapper)
% SEEALSO: vargplvmOptimise
% COPYRIGHT: Andreas C. Damianou, 2012
% VARGPLVM


% model, pruneModel, saveModel, {initVardistIters, itNo}, display
% OR
% model, pruneModel, globalOpt % Take all values from globalOpt struct

pruneModel = true;
saveModel = true;
display = true;

globalOpt = model.globalOpt;

if nargin > 2
    pruneModel = varargin{1};
    if length(varargin) > 1 && isstruct(varargin{2})
        % Take all values from given struct (globalOpt)
        globalOpt = varargin{2};
        saveModel = globalOpt.saveModel;
        display = globalOpt.optimisationDisplay;
    else
        if length(varargin) > 1
            saveModel = varargin{2};
        end
        
        if length(varargin) > 2
            globalOpt.initVardistIters = varargin{3}{1};
            globalOpt.itNo = varargin{3}{2};
        end
        
        if length(varargin) > 3
            display = varargin{4};
        end
    end
end


if ~isfield(model, 'iters')
    model.iters=0;
end
if ~isfield(model, 'initVardistIters')
    model.initVardistIters = 0;
end


i=1;
while ~isempty(globalOpt.initVardistIters(i:end)) || ~isempty(globalOpt.itNo(i:end))
    % do not learn beta for few iterations for intitilization
    if  ~isempty(globalOpt.initVardistIters(i:end)) && globalOpt.initVardistIters(i)
        model.initVardist = 1; model.learnSigmaf = 0;
        fprintf(1,'# Intitiliazing the variational distribution...\n');
        model = vargplvmOptimise(model, display, globalOpt.initVardistIters); % Default: 20
        model.initVardistIters = model.initVardistIters + globalOpt.initVardistIters(i);
        if saveModel
            fprintf('# Saving model after optimising var.distr for %d iterations...\n\n', globalOpt.initVardistIters(i))
            if pruneModel
                modelPruned = vargplvmPruneModel(model);
                vargplvmWriteResult(modelPruned, modelPruned.type, globalOpt.dataSetName, globalOpt.experimentNo);
            else
                vargplvmWriteResult(model, model.type, globalOpt.dataSetName, globalOpt.experimentNo);
            end
        end
    end

    vargplvmShowScales(model, false);
    
    % Optimise the model.
    model.learnBeta = 1;
    model.initVardist = 0;
    model.learnSigmaf = 1;

    if  ~isempty(globalOpt.itNo(i:end)) && globalOpt.itNo(i)
        iters = globalOpt.itNo(i); % Default: 1000
        fprintf(1,'# Optimising the model for %d iterations (session %d)...\n',iters,i);
        model = vargplvmOptimise(model, display, iters);
        if isfield(model,'mOrig')
            fprintf('# SNR = %.6f\n',(1/model.beta)/var(model.mOrig(:)));
        else
            fprintf('# SNR = %.6f\n',(1/model.beta)/var(model.m(:)));
        end
        model.iters = model.iters + iters;
        % Save the results.
        if saveModel
            fprintf(1,'# Saving model after doing %d iterations\n\n',iters)
            if pruneModel
                modelPruned = vargplvmPruneModel(model);
                vargplvmWriteResult(modelPruned, modelPruned.type, globalOpt.dataSetName, globalOpt.experimentNo);
            else
                vargplvmWriteResult(model, model.type, globalOpt.dataSetName, globalOpt.experimentNo);
            end
        end
    end
    i = i+1;
end

