function model = svargplvmModelCreate(arg1, globalOpt, options, optionsDyn)

% SVARGPLVMMODELCREATE Creates a shared VARGPLVM model from a set of VARGPLVM models
% FORMAT
% DESC Creates a shared VARGPLVM model from a set of VARGPLVM models
% ARG model: a cell structure of VARGPLVM models
% RETURN model : the svargplvm model created
%
% SEEALSO : vargplvmCreate
%
% COPYRIGHT : Andreas C. Damianou, 2011

% VARGPLVM

% Check if the first argument is a pre-created model or the training data
% to be used for model creation
if isstruct(arg1{1})
    model = arg1;
else
    Ytr = arg1;
end

if nargin > 1
    numModels = length(Ytr);
    
    %----- Balance modality dimensions
    % Check svargplvn_init for balanceModalityDim entry to see
    % a description of what is going on here.
    if ~iscell(globalOpt.balanceModalityDim)
        tmp = repmat(globalOpt.balanceModalityDim, 1, numModels);
        globalOpt.balanceModalityDim = mat2cell(tmp,1,ones(1,size(tmp,2)));
    end
    modalityMapping = {};
    if sum(cell2mat(globalOpt.balanceModalityDim)) > 0
        maxDval = -Inf;
        
        % Find maximum dimensionality
        for i=1:numModels
            if size(Ytr{i},2) > maxDval
                maxDval = size(Ytr{i},2);
            end
        end
        
        % Do the mapping
        for i=1:numModels
            curD = size(Ytr{i},2);
            if curD < maxDval && globalOpt.balanceModalityDim{i}
                fprintf('# Mapping modality %d from %d to %d dimensions!\n', i, curD, maxDval)
                % Fix seed for the random mapping so that it's reproducible
                curSeed = rng;
                rng(123);
                modalityMapping{i} = rand(curD, maxDval);
                rng(curSeed);
                %--
                Ytmp = zeros(size(Ytr{i},1), maxDval);
                for n = 1:size(Ytmp,1)
                    Ytmp(n,:) = Ytr{i}(n,:)*modalityMapping{i};
                end
                Ytr{i} = Ytmp;
            end
        end
    end


    %-------------- INIT LATENT SPACE ---%
    [X_init, m] = svargplvmInitLatentSpace2(Ytr, globalOpt, options);
    
    if nargin > 3 && ~isempty(globalOpt.dynamicsConstrainType) && ~isempty(optionsDyn)
        dynUsed = 1;
    else
        dynUsed = 0;
    end
    
    
    latentDim = size(X_init,2);
    %-- Create the sub-models: Assume that for each dataset we have one model.
    % This can be changed later, as long as we find a reasonable way to
    % initialise the latent spaces.
    M = length(Ytr);
    for i=1:M
        d{i} = size(Ytr{i},2);
        %---- Here put some code to assign X to the global common X which must
        % be created by doing pca in the concatenation of Y's...After this
        % point, model{i}.X will be the same for all i's. TODO...
        % fprintf(1,'# Creating the model...\n');
        options{i}.initX = X_init;
        model{i} = vargplvmCreate(latentDim, d{i}, Ytr{i}, options{i});
        
        model{i}.X = X_init; %%%%%%%
        model{i} = vargplvmParamInit(model{i}, m{i}, model{i}.X);
        model{i}.X = X_init; %%%%%%%
        var_mi = var(m{i}(:));
        m{i} = []; % Save some memory
        
        if isfield(globalOpt, 'inputScales') && ~isempty(globalOpt.inputScales)
            inpScales = globalOpt.inputScales;
        else
            inpScales = globalOpt.inverseWidthMult./(((max(model{i}.X)-min(model{i}.X))).^2); % Default 5
            %inpScales(:) = max(inpScales); % Optional!!!!!
        end
        
        model{i}.kern.comp{1}.inputScales = inpScales;
        if ~iscell(model{i}.kern)
            model{i}.kern.inputScales = model{i}.kern.comp{1}.inputScales;
        end
        
        
        params = vargplvmExtractParam(model{i});
        model{i} = vargplvmExpandParam(model{i}, params);
        model{i}.vardist.covars = 0.5*ones(size(model{i}.vardist.covars)) + 0.001*randn(size(model{i}.vardist.covars));
        
        
        %-------- Add dynamics to the model -----
        if  dynUsed
            optionsDyn.initX = X_init;
            model{i} = svargplvmAddDynamics(model{i}, globalOpt, optionsDyn);
        end
        
        
        % model{i}.beta=1/(0.01*var_mi);
        % NEW!!!!!
        if var_mi < 1e-8
            warning(['Variance in model ' num2str(i) ' was too small. Setting beta to 1e+7'])
            model{i}.beta = 1e+7;
        else
            model{i}.beta = 1/((1/globalOpt.initSNR * var_mi));
        end
        %        model{i}.beta = 1/((1/globalOpt.initSNR * var(model{i}.m(:)))); %%%%%%%%
        
        %prunedModelInit{i} = vargplvmPruneModel(model{i});
        %disp(model{i}.vardist.covars)
        
        Ytr{i} = []; % Save some memory
    end
end



%---




for i=1:length(model)
    modelNew.comp{i} = model{i};
    modelNew.comp{i}.id = i;
    if isfield(modelNew.comp{i}, 'dynamics') & ~isempty(modelNew.comp{i}.dynamics)
        modelNew.comp{i}.nPrivateParams = modelNew.comp{i}.nParams - modelNew.comp{i}.dynamics.nParams;
    else
        modelNew.comp{i}.nPrivateParams = modelNew.comp{i}.nParams - modelNew.comp{i}.vardist.nParams;
    end
end


% Now set the common fields, they should be the same for any comp{i} at
% this point. That is, the variational distribution and N and q. Also, if
% there are dynamics and since mu = Kt*mubar (and I want mu's to be the
% same) all Kt{i} must be the same. So, in case there are dynamics, the
% whole dynamics structure should be the same.
modelNew.N = modelNew.comp{1}.N;
modelNew.q = modelNew.comp{1}.q;
modelNew.vardist = modelNew.comp{1}.vardist;

% Assume that either all sub-models have dynamics, or none.
if isfield(modelNew.comp{1}, 'dynamics') & ~isempty(modelNew.comp{1}.dynamics)
    modelNew.dynamics = modelNew.comp{1}.dynamics;
end
modelNew.X = modelNew.comp{1}.X;

modelNew.numModels = length(model);

model = modelNew;

% The variational means etc are initialised with adding some random noise.
% If we want the shared parameters to be exactly same for all models, we
% can impose it here.
for i=1:model.numModels
    model.comp{i}.X = model.X;
    model.comp{i}.vardist = model.vardist;
    if isfield(model, 'dynamics') & ~isempty(model.dynamics)
        model.comp{i}.dynamics = model.dynamics;
    end
end

model.type = 'svargplvm';


%--- NEW
% The indices of the parameter vector for the dynamics kernel that exist in
% the globalOpt.fixedKernVarianceIndices field, will stay the same during
% opimisation (i.e their gradients will be forced to be zero)
if nargin>1 && isfield(globalOpt, 'fixedKernVarianceIndices') && ~isempty(globalOpt.fixedKernVarianceIndices)
    if isfield(model, 'dynamics') && ~isempty(model.dynamics)
        model.dynamics.fixedKernVariance = globalOpt.fixedKernVarianceIndices;
        for i=1:length(model.comp)
            model.comp{i}.dynamics.fixedKernVariance = globalOpt.fixedKernVarianceIndices;
        end
    end
end
%----

% % The following field holds the indices that are shared for all sub-models.
% % In the svargplvmExtractParam function, these indices must be ignored
% % because they are the same for all models.
% if isfield(modelNew.comp{1}, 'dynamics') & ~isempty(modelNew.comp{1}.dynamics)
%     model.sharedParams = 1:model.dynamics.nParams;
% else
%     model.sharedParams = 1:model.vardist.nParams;
% end

% Make sure that model.comp{i}.KLweight and dynamics are not used together
if isfield(model.comp{1}, 'dynamics') && ~isempty(model.comp{1}.dynamics)
    for i=1:model.numModels
        if isfield(model.comp{i}, 'KLweight')
            assert(model.comp{1}.KLweight == 0.5)
        end
    end
end


if nargin > 1
    model.optimiser = globalOpt.optimiser;
    
    
    % WARNING: This is not a good thing to do!! See the balanceModalityDim flag
    % for a better solution!
    % Adjust the "influence" of each of the partial likelihood bounds according
    % to their dimensionality. Specifically:
    % Bound is: F = F_1 + F_2 + ... + F_M - KL
    % a_i = 1/model.comp{i}.d = 1 / di
    % First we balance the bound with the dimensionality by multiplying each
    % F_i with its corresponding a_i
    % But now the KL is unequally balanced, since overall we might have given
    % more importance to the F_i terms (in terms of gradient steps, since all
    % multipliers are just scales of the corresponding partial grad). So, we
    % then multiply everything with d1*d2*d3*... / (d2*d3+d1*d3+d1*d2+..). In
    % this way we make sure that overall sum(F_i) coefficients is 1, as if we
    % didn't weight anything. Notice that the weights a_i are given by
    % 1-(a_i/2) (the KL is not actually affected as it is
    % computed in the svargplvm and not the vargplvm framework) because in the
    % vargplvm framework the weights are automatically multiplied by 2.
    if globalOpt.balanceModalities
        % I know this is a stupid way of computing it, but I want to go for
        % lunch.
        fprintf('# Balancing modalities...\n');
        dimProd = 1;
        for i=1:model.numModels
            dimProd = dimProd * log(model.comp{i}.d);
        end
        for i=1:model.numModels
            dimPair(i) = dimProd/log(model.comp{i}.d);
        end
        dimPairSum = sum(dimPair);
        for i=1:model.numModels
            a_i = (dimProd / dimPairSum)/log(model.comp{i}.d);
            model.comp{i}.KLweight = 1-(a_i/2);
            assert(model.comp{i}.KLweight <= 1 && model.comp{i}.KLweight>=0)
        end
    end
    
    %{
    if ~isempty(modalityMapping)
        % This is too costly, don't store...
        model.modalityMapping = modalityMapping;
        %Instead, the random seed to generate it is
        % fixed in this function, so modalityMapping{i} can be reproduced by:
        %curSeed = rng;
        %rng(1);
        %modalityMapping{i} = rand(model.comp{i}.d, maxDval);
        %rng(curSeed);
    end
    %}
end