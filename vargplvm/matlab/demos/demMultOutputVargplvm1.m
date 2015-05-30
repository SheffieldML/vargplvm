% DEMMULTOUTPUTVARGPLVM1 Run a shared var-GPLVM with D sub-models which is equivalent to a vargplvm with multiple outputs.
% DESC Run a shared var-GPLVM with D sub-models which is equivalent to a vargplvm with multiple outputs.
%
% COPYRIGHT: Andreas C. Damianou, Carl Henrik Ek, 2011
% SEEALSO: demSharedVargplvm1
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
if ~exist('mappingKern')   ,  mappingKern = {'rbfard2', 'white'}; end
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
if ~exist('initLatent'),     initLatent ='ppcaConcatenate';   end % That's for the whole latent space init.
if ~exist('dataToKeep'), dataToKeep = -1 ;end
if ~exist('toyDataCreate'), toyDataCreate = 'vargplvm'; end
if ~exist('doPredictions'), doPredictions = 0; end
if ~exist('dataType'), dataType = []; end

%%


%load humanPose;
Ytmp = svargplvmLoadData('humanPose');
Z=Ytmp{3};
clear 'Ytmp';
%ind2 = floor(1:2:size(Y,1));
%%%__silhouettes are whole images
if exist('imageSilhouette') && imageSilhouette
    load sil_images
end
% Remove the 'drunk' sequence
Z1 = Z(1:100,:);
Z2 = Z(337:end,:);
Z = [Z1; Z2];
clear('Z1','Z2','Y');

% Subsample
ind2 = floor(1:2:size(Z,1));
%tmp=setdiff(1:size(Z,1),ind2);
Z = Z(ind2,:);
if exist('dataToKeep') & dataToKeep ~= -1
    Z = Z(1:dataToKeep,:);
end
dataSetNames = {};
for d=7:size(Z,2)
    Y{d-6} = Z(:,d);
    dataSetNames={dataSetNames{:}, ['data' num2str(d-6)]};
end


clear d
%xyzankurAnim(Z, 3);

numberOfDatasets = length(Y);
%initLatent = 'ppcaConcatenate';
latentDimPerModel=1;
numSubModels = numberOfDatasets;
numSharedDims = 0;

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
    options{i}.optimiser = 'scg';
    
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
%!!!! Doesn't work if Y{i}'s have different sizes!!
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
    X_init = initLatent;
elseif strcmp(initLatent, 'ncca')
    %-- Learn Initialisation through NCCA ( FOR TWO DATASETS only) %%%!!!
    if size(Ytr) ~= 2
        error('ncca initialization only when there are two datasets!');
    end
    [Xsy Xsz Xy Xz] = nccaEmbed(Ytr{1},Ytr{2},uint8([7 7]),uint8(1),uint8([2 2]),true);
    Xs = (1/2).*(Xsy+Xsz);
    X_init = [Xy Xs Xz]; % sizes: 2,1,2
    X_init = (X_init-repmat(mean(X_init),size(X_init,1),1))./repmat(std(X_init),size(X_init,1),1);
elseif strcmp(initLatent,'ppca')
    %-- Learn initialisation through PCA: Perform mappings from Y_i to X_i
    % and concatenate X_i's to augment the X's dimensionality as more
    % datasets are added.
    %     X_init = [];
    %     for i=1:numberOfDatasets
    %         initFunc = str2func([initX 'Embed']);
    %         X_init_cur = initFunc(m{i}, latentDimPerModel);
    %         X_init = [X_init X_init_cur];
    %     end
    X_init{1} = ppcaEmbed(m{1},7);
    X_init{2} = ppcaEmbed(m{2},3);
    X_init = [X_init{1} X_init{2}];
elseif strcmp(initLatent,'ppcaConcatenate')
    initFunc = str2func([initX 'Embed']);
    X_init = initFunc(mAll, latentDimPerModel * numSubModels + numSharedDims);
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
elseif strcmp(initLatent, 'outputs')
    X_init = mAll;
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
    model{i}.X = X_init;
    model{i} = vargplvmParamInit(model{i}, model{i}.m, model{i}.X);
    model{i}.X = X_init;
    
    inpScales = invWidthMult./(((max(model{i}.X)-min(model{i}.X))).^2); % Default 5
    
    %inpScales(:) = max(inpScales); % Optional!!!!!
    if strcmp(model{i}.kern.type, 'cmpnd') && (strcmp(model{i}.kern.comp{1}.type, 'rbfard2') || strcmp(model{i}.kern.comp{1}.type, 'linard2'))
        model{i}.kern.comp{1}.inputScales = inpScales;
    elseif strcmp(model{i}.kern.type, 'rbfardjit')
        model{i}.kern.inputScales = inpScales;
    end

    params = vargplvmExtractParam(model{i});
    model{i} = vargplvmExpandParam(model{i}, params);
    
    
    model{i}.vardist.covars = 0.5*ones(size(model{i}.vardist.covars)) + 0.001*randn(size(model{i}.vardist.covars));
    
    %-------- Add dynamics to the model -----
 %{
    if dynUsed
        fprintf(1,'# Adding dynamics to the model...\n');
        optionsDyn{i}.type = 'vargpTime';
        optionsDyn{i}.t=timeStampsTraining{i};
        optionsDyn{i}.inverseWidth=invWidthMultDyn; % Default: 100
        optionsDyn{i}.initX = X_init; % initX; % That's for the dynamics
        
        kern = kernCreate(t{i}, dynamicKern); % Default: {'rbf','white','bias'}
        
        %-- Initialize each element of the compound kernel
        % ATTENTION: For the gradients we assume that the base kernel (rbf,
        % matern etc) must be the FIRST one and if a second base kernel
        % (not white or bias) exist must be the LAST one!!!!!!!!!!!!!!
        if isfield(kern,'comp')
            fprintf('# Dynamics Kernel initialization: \n')
            kernFound = 0;
            for k=1:numel(kern.comp)
                type = kern.comp{i}.type;
                if strcmp(type, 'rbfperiodic') || strcmp(type, 'rbfperiodic2')
                    if exist('periodicPeriod')
                        kern.comp{k}.period = periodicPeriod;
                        kern.comp{k}.factor = 2*pi/periodicPeriod;
                    end
                    fprintf(1,'\t # periodic period: %d\n',kern.comp{k}.period);
                elseif strcmp(type, 'whitefixed')
                    if ~exist('whiteVar')
                        whiteVar = 1e-6;
                    end
                    kern.comp{k}.variance = whiteVar;
                    fprintf(1,'\t # fixedwhite variance: %d\n',whiteVar);
                elseif strcmp(type, 'white')
                    if ~exist('whiteVar')
                        %     whiteVar = 1e-4; % Some models have been trained
                        %     with this!!
                        whiteVar = 0.1;
                    end
                    fprintf(1,'\t # white variance: %d\n',whiteVar);
                    kern.comp{k}.variance = whiteVar; % Usual values: 1e-1, 1e-3
                elseif strcmp(type, 'bias')
                    if exist('biasVar')
                        kern.comp{k}.bias = biasVar;
                        fprintf('\t # bias variance: %d \n', biasVar);
                    end
                end
                % The following is related to the expected number of
                % zero-crossings.(larger inv.width numerator, rougher func)
                if strcmp(type,'rbfperiodic') || strcmp(type,'rbfperiodic2') || strcmp(type,'rbf') || strcmp(type,'matern32')
                    kern.comp{k}.inverseWidth = optionsDyn{i}.inverseWidth./(((max(t{i})-min(t{i}))).^2);
                    kern.comp{k}.variance = 1;
                    % This is a bit hacky: if this is the second time an
                    % rbf, or rbfperiodic or... kernel is found, then the
                    % second variance can be initialised to be smaller
                    if kernFound
                        if ~exist('secondVarInit')
                            kern.comp{k}.variance = 0.006;
                        else
                            kern.comp{k}.variance = secondVarInit;
                        end
                        fprintf('\t # Second variance initialized to %d \n',kern.comp{k}.variance);
                        
                    end
                    kernFound = k;
                end
            end
        end
        
        optionsDyn{i}.kern = kern;
        optionsDyn{i}.vardistCovars = vardistCovarsMult; % 0.23 gives true vardist.covars around 0.5 (DEFAULT: 0.23) for the ocean dataset
        
        % Fill in with default values whatever is not already set
        optionsDyn{i} = vargplvmOptionsDyn(optionsDyn{i});
        model{i} = vargplvmAddDynamics(model{i}, 'vargpTime', optionsDyn{i}, optionsDyn{i}.t, 0, 0,optionsDyn{i}.seq);
        
        fprintf(1,'# Further calibration of the initial values...\n');
        model{i} = vargplvmInitDynamics(model{i},optionsDyn{i});
        
        %  to also not learn the last kernel's variance
        if numel(kern.comp) > 1 && exist('learnSecondVariance') && ~learnSecondVariance
            fprintf(1,'# The variance for %s in the dynamics is not learned!\n',kern.comp{end}.type)
            model{i}.dynamics.learnSecondVariance = 0;
            model{i}.dynamics.kern.comp{end}.inverseWidth = model{i}.dynamics.kern.comp{1}.inverseWidth/10; %%% TEMP
        end
    end
%}
    
    model{i}.beta=1/(0.01*var(model{i}.m(:))); % Initial SNR is set to 100
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
dataType = 'humanPose';
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
profile on %%%%%%

display = 1;
%%%% Optimisation
% do not learn beta and sigma_f for few iterations for intitialization
if initVardistIters ~=0
    model.initVardist = 1;
    model = svargplvmPropagateField(model,'initVardist', model.initVardist);
    fprintf(1,'# Intitiliazing the variational distribution...\n');
    model = svargplvmOptimise(model, display, initVardistIters); % Default: 20
    %fprintf(1,'1/b = %.4d\n',1/model.beta);
    
    modelInitVardist = model;
    model.initVardistIters=initVardistIters;
end

model.initVardist = 0;
model = svargplvmPropagateField(model,'initVardist', model.initVardist);

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

profile off %%%%%%%

save(fileToSave, 'model', 'prunedModelInit');
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