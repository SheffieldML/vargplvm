function model = svargplvmOptimiseModelNoDisplay(model, varargin)
% SVARGPLVMOPTIMISEMODELNODISPLAY See svargplvmOptimiseModel
%
% VARGPLVM

% modelType can be 'bc' or 'discr'


% model, pruneModel, saveModel, {initVardistIters, itNo}

pruneModel = false;
saveModel = true;

globalOpt = model.globalOpt;


if nargin > 2
    pruneModel = varargin{1};
    if length(varargin) > 1
        saveModel = varargin{2};
    end

    if length(varargin) > 2
        globalOpt.initVardistIters = varargin{3}{1};
        globalOpt.itNo = varargin{3}{2};
    end
end



model.iters=0;
model.initVardistIters = 0;
display = 0;


i=1;
while ~isempty(globalOpt.initVardistIters(i:end)) || ~isempty(globalOpt.itNo(i:end))
    % do not learn beta for few iterations for intitilization
    if  ~isempty(globalOpt.initVardistIters(i:end)) && globalOpt.initVardistIters(i)
        model.initVardist = 1; model.learnSigmaf = 0;
        model = svargplvmPropagateField(model,'initVardist', model.initVardist);
        model = svargplvmPropagateField(model,'learnSigmaf', model.learnSigmaf);
        fprintf(1,'# Intitiliazing the variational distribution...\n');
        model = svargplvmOptimise(model, display, globalOpt.initVardistIters); % Default: 20
        model.initVardistIters = model.initVardistIters + globalOpt.initVardistIters(i);
        if saveModel
            fprintf('# Saving model after optimising beta for %d iterations...\n\n', globalOpt.initVardistIters(i))
            if pruneModel
                modelPruned = svargplvmPruneModel(model);
                vargplvmWriteResult(modelPruned, modelPruned.type, globalOpt.dataSetName, globalOpt.experimentNo);
            else
                vargplvmWriteResult(model, model.type, globalOpt.dataSetName, globalOpt.experimentNo);
            end
        end
    end

    % Optimise the model.
    model.learnBeta = 1;
    model.initVardist = 0;
    model.learnSigmaf = 1;

    model.date = date;
    if  ~isempty(globalOpt.itNo(i:end)) && globalOpt.itNo(i)
        model = svargplvmPropagateField(model,'initVardist', model.initVardist);
        model = svargplvmPropagateField(model,'learnBeta', model.learnBeta);
        model = svargplvmPropagateField(model,'learnSigmaf', model.learnSigmaf);

        iters = globalOpt.itNo(i); % Default: 1000
        fprintf(1,'# Optimising the model for %d iterations (session %d)...\n',iters,i);
        model = svargplvmOptimise(model, display, iters);
        for j=1:length(model.comp)
            fprintf('# SNR%d = %.6f\n',j,(1/model.comp{j}.beta)/var(model.comp{j}.m(:)));
        end
        SNR = svargplvmSNR(model);
        svargplvmCheckSNR(SNR);
        
        model.iters = model.iters + iters;
        % Save the results.
        if saveModel
            fprintf(1,'# Saving model after doing %d iterations\n\n',iters)
            if pruneModel
                modelPruned = svargplvmPruneModel(model);
                vargplvmWriteResult(modelPruned, modelPruned.type, globalOpt.dataSetName, globalOpt.experimentNo);
            else
                vargplvmWriteResult(model, model.type, globalOpt.dataSetName, globalOpt.experimentNo);
            end
        end
    end
    i = i+1;
end

