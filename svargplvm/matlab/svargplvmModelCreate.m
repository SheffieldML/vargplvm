function model = svargplvmModelCreate(arg1, globalOpt, options, optionsDyn)

% SHEFFIELDMLMODELCREATE Creates a shared VARGPLVM model from a set of VARGPLVM models
% FORMAT
% DESC Creates a shared VARGPLVM model from a set of VARGPLVM models
% ARG m odel: a cell structure of VARGPLVM models
% RETURN model : the svargplvm model created
%
% SEEALSO : vargplvmCreate
%
% COPYRIGHT : Andreas C. Damianou, 2011

% SHEFFIELDML

% Check if the first argument is a pre-created model or the training data
% to be used for model creation
if isstruct(arg1{1})
    model = arg1;
else
    Ytr = arg1;
end

if nargin > 1
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
    for i=1:length(Ytr)
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


        % model{i}.beta=1/(0.01*var(m{i}(:)));
        % NEW!!!!!
        if var(m{i}(:)) < 1e-8
            warning(['Variance in model ' num2str(i) ' was too small. Setting beta to 1e+7'])
            model{i}.beta = 1e+7;
        else
            model{i}.beta = 1/((1/globalOpt.initSNR * var(m{i}(:))));
        end
%        model{i}.beta = 1/((1/globalOpt.initSNR * var(model{i}.m(:)))); %%%%%%%%
        
        %prunedModelInit{i} = vargplvmPruneModel(model{i});
        %disp(model{i}.vardist.covars)
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
