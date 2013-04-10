% DEMHIERARCHICAL A simple hierarchical demo, where we run a
% multiple-output svargplvm on the latent space found by running a normal
% svargplvm on the Yale faces. The multiple-output svargplvm puts one
% shared sub-model for each latent dimension.
% TODO: Visualisation

% VARGPLVM



% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

% Define constants (in a manner that allows other scripts to parametrize
% this one).
if ~exist('experimentNo') ,  experimentNo = 4044;      end
if ~exist('itNo')         ,  itNo = [500 500];              end     % Default: 2000
if ~exist('indPoints')    ,  indPoints = 80;          end     % Default: 49
if ~exist('latentDimPerModel')    ,  latentDimPerModel = 3;          end
% Set to 1 to use dynamics or to 0 to use the standard var-GPLVM
if ~exist('dynUsed')      ,  dynUsed = 0;             end
if ~exist('initVardistIters'), initVardistIters = 100;      end
if ~exist('dynamicKern')   ,     dynamicKern = {'rbf', 'white', 'bias'}; end
if ~exist('mappingKern')   ,  mappingKern = 'rbfardjit'; end % {'rbfard2', 'white'}; end
% if ~exist('mappingKern'),  mappingKern = {'rbfard2', 'bias', 'white'}; end
% Set to 1 to tie the inducing points with the latent vars. X
if ~exist('fixInd')        ,     fixInd = 0;    end
% 0.1 gives around 0.5 init.covars. 1.3 biases towards 0.
if ~exist('vardistCovarsMult'),  vardistCovarsMult=1.3;                  end
% Set to empty value {} to work with toy data
if ~exist('dataSetNames')    ,    dataSetNames = {'USecon','NYSE2Small'};    end
if ~exist('invWidthMultDyn'),    invWidthMultDyn = 100;                     end
if ~exist('invWidthMult'),       invWidthMult = 5;                     end
if ~exist('initX'),     initX ='ppca';   end % That's for the dynamics initialisation
if ~exist('initLatent'),     initLatent = 'ppcaConcatenate';   end % That's for the whole latent space init.
if ~exist('dataToKeep'), dataToKeep = -1 ;end
if ~exist('toyDataCreate'), toyDataCreate = 'vargplvm'; end
if ~exist('doPredictions'), doPredictions = 0; end
if ~exist('dataType'), dataType = []; end
if ~exist('latentDim'), latentDim = 14; end

%%
dataType = 'Yale';


load Xexp25.mat
Yconcatenated = X;


for d=1:size(Yconcatenated,2)
    Y{d} = Yconcatenated(:,d);
end

clear('d','X');

numberOfDatasets = length(Y);
latentDimPerModel=1;
numSubModels = numberOfDatasets;

%--


%%

learnScales = 0;
mAll=[];

%-- Load datasets
for i=1:numberOfDatasets
    Y_cur = Y{i};
    dims{i} = size(Y_cur,2);
    N{i} = size(Y_cur,1);
    indTr = 1:N{i};
    indTs = setdiff(size(Y_cur,1), indTr);
    Ytr{i} = Y_cur(indTr,:); %Yts = Y(indTs,:);
    d{i} = size(Ytr{i}, 2);
end

%%%%%% TEMP: N{i}'s must be the same!!
for i=2:numberOfDatasets
    if N{i} ~= N{i-1}
        error('The number of observations in each dataset must be the same!');
    end
end
%%%%


%-- Options for the models
for i=1:numberOfDatasets
    % Set up models
    options{i} = vargplvmOptions('dtcvar');
    options{i}.kern = mappingKern; %{'rbfard2', 'bias', 'white'};
    %indPoints = 80; %%%%%
    options{i}.numActive = indPoints;
    options{i}.optimiser = 'scg2';
    
    % !!!!! Be careful to use the same type of scaling and bias for all
    % models!!!
    
    % scale = std(Ytr);
    % scale(find(scale==0)) = 1;
    %options.scaleVal = mean(std(Ytr));
    options{i}.scaleVal = sqrt(var(Ytr{i}(:)));
    if fixInd
        options{i}.fixInducing=1;
        options{i}.fixIndices=1:size(Ytr{i},1);
    end
end


%-- Create the normalised version of the datasets and concatenate
for i=1:numberOfDatasets
    % Compute m, the normalised version of Ytr (to be used for
    % initialisation of X)
    bias = mean(Ytr{i});
    scale = ones(1, d{i});
    
    if(isfield(options{i},'scale2var1'))
        if(options{i}.scale2var1)
            scale = std(Ytr{i});
            scale(find(scale==0)) = 1;
            if(learnScales)
                warning('Both learn scales and scale2var1 set for GP');
            end
            if(isfield(options{i}, 'scaleVal'))
                warning('Both scale2var1 and scaleVal set for GP');
            end
        end
    end
    if(isfield(options{i}, 'scaleVal'))
        scale = repmat(options{i}.scaleVal, 1, d{i});
    end
    
    % Remove bias and apply scale.
    m{i} = Ytr{i};
    for j = 1:d{i}
        m{i}(:, j) = m{i}(:, j) - bias(j);
        if scale(j)
            m{i}(:, j) = m{i}(:, j)/scale(j);
        end
    end
    
    mAll = [mAll m{i}]; % Concatenation (doesn't work if different sizes)
end

% Clear some variables
clear('Y','Y_cur','bias','scale','ind2');


% %-- Create shared X:
% initFunc = str2func([initX 'Embed']);
% X = initFunc(mAll, latentDim);
if ~isstr(initLatent)
     % here, consider an initialisation of the shared and two private parts
     % separately
    X_init = initLatent;
elseif strcmp(initLatent,'same')
    if latentDim ~= size(Yconcatenated,2)
        X_init = ppcaEmbed(Yconcatenated, latentDim);
    else
        X_init = Yconcatenated;
    end
elseif strcmp(initLatent,'ppcaConcatenate')
    initFunc = str2func([initX 'Embed']);
    X_init = initFunc(mAll, latentDim);
elseif strcmp(initLatent, 'pca')
    % We like the numer of latent dims to be 2+numSharedDims, ideally 3. With
    % vargplvm we set Q=6 and expect the other 3 or 4 to be switched off.
    [U,V] = pca(mAll,latentDim);
    X_init = mAll*V;
elseif strcmp(initLatent, 'pca2')
    % We like the numer of latent dims to be 2+numSharedDims, ideally 3. With
    % vargplvm we set Q=6 and expect the other 3 or 4 to be switched off.
    [U,V] = pca(m{1},latentDim);
    X_init{1} = m{1}*V;
    [U,V] = pca(m{2},latentDim);
    X_init{2} = m{2}*V;
    X_init = [X_init{1} X_init{2}];
elseif strcmp(initLatent, 'pca3')
    clear mAll
    % We like the number of latent dims to be 2+numSharedDims, ideally 3. With
    % vargplvm we set Q=6 and expect the other 3 or 4 to be switched off.
    try
        [U,V] = pca(m{1},7);
        X_init{1} = m{1}*V;
        m{1}=[];%%%
        [U,V] = pca(m{2},3);
        X_init{2} = m{2}*V;
        X_init = [X_init{1} X_init{2}];
    catch e
        if strcmp(e.identifier, 'MATLAB:nomem')
            fprintf('# !!! Warning: Not enough memory to initialise with PCA! Initialising with %s instead...\n',initX);
        end
        initFunc = str2func([initX 'Embed']);
        X_init{1} = initFunc(m{1}, 7);
        X_init{2} = initFunc(m{2},3);
        X_init = [X_init{1} X_init{2}];
    end
end
latentDim = size(X_init,2);

% Free up some memory
clear('m')



%-- Create the sub-models: Assume that for each dataset we have one model.
% This can be changed later, as long as we find a reasonable way to
% initialise the latent spaces.
for i=1:numberOfDatasets
    %---- Here put some code to assign X to the global common X which must
    % be created by doing pca in the concatenation of Y's...After this
    % point, model{i}.X will be the same for all i's. TODO...
    fprintf(1,'# Creating the model...\n');
    options{i}.initX = X_init;
    model{i} = vargplvmCreate(latentDim, d{i}, Ytr{i}, options{i});
    model{i}.X = X_init; %%%%%%%
    model{i} = vargplvmParamInit(model{i}, model{i}.m, model{i}.X);
    model{i}.X = X_init; %%%%%%%
    
    inpScales = invWidthMult./(((max(model{i}.X)-min(model{i}.X))).^2); % Default 5
    
    %inpScales(:) = max(inpScales); % Optional!!!!!
    
    model{i}.kern.comp{1}.inputScales = inpScales;
    
    if strcmp(model{i}.kern.type, 'rbfardjit')
        model{i}.kern.inputScales = model{i}.kern.comp{1}.inputScales;
    end
    params = vargplvmExtractParam(model{i});
    model{i} = vargplvmExpandParam(model{i}, params);
    
    
    model{i}.vardist.covars = 0.5*ones(size(model{i}.vardist.covars)) + 0.001*randn(size(model{i}.vardist.covars));
 
    
    model{i}.beta=1/(0.01*var(model{i}.m(:)));
    prunedModelInit{i} = vargplvmPruneModel(model{i});
    %disp(model{i}.vardist.covars)
end

% if experimentNo == -1
%     experimentNo = globalExperimentNumber('sharedVargplvm', 'dataType');
% end


%modelInit = model;%%%TEMP

% TODO:
%--  Unify models into a structure
model = svargplvmModelCreate(model);
model.dataSetNames = dataSetNames;
model.initLatent = initLatent;
model.experimentNo = experimentNo;
model.dataType = dataType;
%%---
capName = dataType;
capName(1) = upper(capName(1));
modelType = model.type;
modelType(1) = upper(modelType(1));
fileToSave = ['dem' capName modelType num2str(experimentNo) '.mat'];
%%---

%-- Define what level of parallelism to use (w.r.t submodels or/and w.r.t
% datapoints).
%if model.numModels > 8
    fprintf('# Parallel computations w.r.t the submodels!\n');
    model.parallel = 1;
    model = svargplvmPropagateField(model,'parallel', 1);
%elseif model.N > 15
%    fprintf('# Parallel computations w.r.t the datapoints!\n');
%    model.vardist.parallel = 1;
%    for i=1:model.numModels
%        model.comp{i}.vardist.parallel = 1;
%    end
%end


%%%
params = svargplvmExtractParam(model);
model = svargplvmExpandParam(model, params);
%    model = svargplvmOptimise(model, 1, 200);

%%
% %%%
%  for i=1:length(model.comp)
%     subplot(4,ceil(numberOfDatasets/4),i)
%     bar(model.comp{i}.kern.inputScales);
%  end
%%%

display = 1;
%%%% Optimisation
% do not learn beta and sigma_f for few iterations for intitialization
if initVardistIters ~=0
    model.initVardist = 1; model.learnSigmaf =0; model.learnBeta = 0;
    model = svargplvmPropagateField(model,'initVardist', model.initVardist);
    model = svargplvmPropagateField(model,'learnSigmaf', model.learnSigmaf);
    model = svargplvmPropagateField(model,'learnBeta', model.learnBeta);
    fprintf(1,'# Intitiliazing the variational distribution...\n');
    model = svargplvmOptimise(model, display, initVardistIters); % Default: 20
    %fprintf(1,'1/b = %.4d\n',1/model.beta);
    
    modelInitVardist = model;
    model.initVardistIters=initVardistIters;
end

model.initVardist = 0; model.learnSigmaf = 1; model.learnBeta = 1;
model = svargplvmPropagateField(model,'initVardist', model.initVardist);
model = svargplvmPropagateField(model,'learnSigmaf', model.learnSigmaf);
model = svargplvmPropagateField(model,'learnBeta', model.learnBeta);

model.iters = 0;



% Optimise the model.
for i=1:length(itNo)
    iters = itNo(i); % default: 2000
    fprintf(1,'\n# Optimising the model for %d iterations (session %d)...\n',iters,i);
    model = svargplvmOptimise(model, display, iters);
    model.iters = model.iters + iters;
    % fprintf(1,'1/b = %.4d\n',1/model.beta);
    %fprintf(1,'# 1/b=%.4f\n var(m)=%.4f\n',1/model.beta, var(model.mOrig(:)));
    
    % Save model
    fprintf(1,'# Saving %s\n',fileToSave);
    %prunedModel = vargplvmPruneModel(model);
    % prunedModelTr = vargplvmPruneModel(modelTr);
    save(fileToSave, 'model', 'prunedModelInit');
end


save(fileToSave, 'model', 'prunedModelInit');

for i=1:length(Ytr)
    subplot(4,4,i)
    bar(model.comp{i}.kern.inputScales)
    fprintf('SNR #%d = %.5f\n', i, var(model.comp{i}.m(:))/model.comp{i}.beta)
end


%{
%prunedModelTr = prunedModel;
%save(fileToSave, 'model', 'prunedModelInit', 'prunedModelTr');

% for i=1:numberOfDatasets
%     figure, bar(prunedModelInit{i}.kern.comp{1}.inputScales);
%     title(['Init scales for dataset ' num2str(i)]);
% end
% for i=1:numberOfDatasets
%     bar(model.comp{i}.kern.comp{1}.inputScales);
%     title(['Final scales for dataset ' num2str(i)]);
%     pause
% end



%%
%--------------------- PREDICTIONS --------------%
% For predictions, we will predict each dimension d separately; namely, we
% will use model.comp{d} to do predictions for Ytest(:,d).

if ~doPredictions
    return
end

% TODO...

%}