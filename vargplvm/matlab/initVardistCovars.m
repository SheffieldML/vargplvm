% This function is not used anywhere...

function best = initVardistCovars
clear
best = -1; bestVal = -1;
vardistCovarsMult = 0;
searchInds =  0.06:0.01:0.2; % Default:  0.08:0.01:0.2 OR 0.06:0.01:0.2
for vardistCovarsMult = searchInds
    disp(vardistCovarsMult)
    try
        model = tempFunc(vardistCovarsMult);
        params = vargplvmExtractParam(model);
        model = vargplvmExpandParam(model, params);
        disp(mean(model.vardist.covars))
        if min(model.vardist.covars) > 0.1 & min(model.vardist.covars) > bestVal
            best = vardistCovarsMult
            bestVal = min(model.vardist.covars)
            fprintf('\n')
        end
    catch e
         fprintf('ooops!!!!\n')
    end
end


%%%%%%%%%% The following has to be replaced with the code that creates and
%%%%%%%%%% initialises the model!!!!!!!!!!!!
function model = tempFunc(vardistCovarsMult)

initial_X = 'separately';
dynUsed=1;
inds = [{1:100} {337:405} {406:471} {1532:1769} {1770:1927}];
initVardistIters = 190;
itNo = [800 700 800 800 900];
indPoints = 100;

%___
% initial_X = 'separately';
% dynUsed=1;
% inds = [{1:100} {101:336} {337:405} {406:471} {1532:1769} {1770:1927}];
% initVardistIters = 190;
% mappingKern = 'rbfardjit';
% itNo = [800 700 800 800 900];
% indPoints = 100;

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

% Define constants (in a manner that allows other scripts to parametrize
% this one).
if ~exist('experimentNo') ,  experimentNo = 404;      end
if ~exist('itNo')         ,  itNo = [500 500];              end     % Default: 2000
if ~exist('indPoints')    ,  indPoints = 100;          end     % Default: 49
if ~exist('latentDimPerModel')    ,  latentDimPerModel = 6;          end
% Set to 1 to use dynamics or to 0 to use the standard var-GPLVM
if ~exist('dynUsed')      ,  dynUsed = 0;             end
if ~exist('initVardistIters'), initVardistIters = 180;      end
if ~exist('dynamicKern')   ,     dynamicKern = {'rbf', 'white', 'bias'}; end
if ~exist('mappingKern')   ,  mappingKern = {'rbfard2', 'white', 'bias'}; end


% 0.1 gives around 0.5 init.covars. 1.3 biases towards 0.
if ~exist('vardistCovarsMult'),  vardistCovarsMult=1.3;                  end
% Set to empty value {} to work with toy data

if ~exist('invWidthMultDyn'),    invWidthMultDyn = 100;                     end
if ~exist('invWidthMult'),       invWidthMult = 5;                     end
if ~exist('initX'),     initX ='ppca';   end % That's for the dynamics initialisation
if ~exist('dataType'), dataType = 'humanPose'; end
if ~exist('enableParallelism'), enableParallelism = 1; end
if ~exist('DgtN'), DgtN = false; end
% Create initial X by doing e.g. ppca in the concatenated model.m's or by
% doing ppca in the model.m's separately and concatenate afterwards?
if ~exist('initial_X'), initial_X = 'separately'; end % Other options: 'together'
% Which indices to use for training, rest for test
if ~exist('indTr'), indTr = -1; end


load sil_images
Yim = Y;
Yim_test = Y_test;

load humanPose;

%ind2 = floor(1:2:size(Y,1));
%%%__silhouettes are whole images. If imageSilhouette is true, then all
%%%pixels will be used (instead of just features extracted).
if exist('imageSilhouette') && imageSilhouette
    fprintf('# Using the actual silhouette''s pixels!\n');
    Y = Yim;
    Y_test = Yim_test;
    %clear Yim
    %clear Yim_test
end


% Subsample
Yall{1} = Y;
Yall{2} = Z;

%---
subSamp = 1;
if exist('inds')
    seqOrig = seq;
    %inds = [{1:100} {337:405} {406:471} {955:1015}];
    seq = []; indsAll = []; prev = 0;
    for i=1:length(inds)
        inds{i} = inds{i}(1:subSamp:end);
        indsAll = [indsAll inds{i}];
        prev = prev + length(inds{i});
        seq = [seq prev];
    end
    Yall{1} = Yall{1}(indsAll,:);
    Yall{2} = Yall{2}(indsAll,:);
    Yim = Yim(indsAll,:);
end

dataSetNames={'silhouette', 'pose'};
% X_init = Xp;
%mappingKern = {'linard2', 'white'};
%mappingKern = {'rbfard2', 'white'};
latentDim = 5; % Anything > 2 and < 10
%xyzankurAnim(Z_test, 3);
clear Y
numberOfDatasets = length(Yall);



%-- Load datasets
for i=1:numberOfDatasets
    Y = Yall{i};
    
    dims{i} = size(Y,2);
    N{i} = size(Y,1);
    if indTr == -1
        indTr = 1:N{i};
    end
    t{i} = linspace(0, 2*pi, size(Yall{i}, 1)+1)'; t{i} = t{i}(1:end-1, 1);
    
    indTs = setdiff(1:size(Y,1), indTr);
    Ytr{i} = Y(indTr,:);
    Yts{i} = Y(indTs,:);
    
    d{i} = size(Ytr{i}, 2);
end
timeStampsTraining = t{1}(indTr,1); %timeStampsTest = t(indTs,1);


for i=2:numberOfDatasets
    if N{i} ~= N{i-1}
        error('The number of observations in each dataset must be the same!');
    end
end


%{
for i=1:size(Ytr{1},1)
    %figure
    %subplot(1,2,2)
    %imagesc(reshape(Ytr{1}(i,:), height, width)), colormap('gray');
    %subplot(1,2,1)
    handle = xyzankurVisualise(Ytr{2}(i,:));
    pause;
    close
    i
end
%}

%%

%-- Options for the models
for i=1:numberOfDatasets
    % Set up models
    options{i} = vargplvmOptions('dtcvar');
    options{i}.kern = mappingKern; %{'rbfard2', 'bias', 'white'};
    options{i}.numActive = indPoints;
    options{i}.optimiser = 'scg2';
    if ~DgtN
        options{i}.enableDgtN = false;
    end
    % !!!!! Be careful to use the same type of scaling and bias for all
    % models!!!
    
    % scale = std(Ytr);
    % scale(find(scale==0)) = 1;
    %options.scaleVal = mean(std(Ytr));
    options{i}.scaleVal = sqrt(var(Ytr{i}(:)));
    if dynUsed
        optionsDyn{i}.type = 'vargpTime';
        optionsDyn{i}.t=timeStampsTraining;
        optionsDyn{i}.inverseWidth=30;
        optionsDyn{i}.seq = seq;
        % Fill in with default values whatever is not already set
        optionsDyn{i} = vargplvmOptionsDyn(optionsDyn{i});
    end
end

%-------------- INIT LATENT SPACE ---%
for i=1:length(Ytr)
    % Compute m, the normalised version of Ytr (to be used for
    % initialisation of X)
    bias = mean(Ytr{i});
    scale = ones(1, d{i});
    
    if(isfield(options{i},'scale2var1'))
        if(options{i}.scale2var1)
            scale = std(Ytr{i});
            scale(find(scale==0)) = 1;
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
    
    %    mAll = [mAll m{i}]; % Concatenation (doesn't work if different sizes)
end





if strcmp(initial_X, 'separately')
    fprintf('# Initialising X by performing ppca in each observed (scaled) dataset separately and then concatenating...\n');
%     X_init{1} = ppcaEmbed(m{1},latentDimPerModel);
%     X_init{2} = ppcaEmbed(m{2},latentDimPerModel);
    X_init{1} = ppcaEmbed(m{1},7);
    X_init{2} = ppcaEmbed(m{2},3);
    X_init = [X_init{1} X_init{2}];
%     
%      [U,V] = pca(m{1},7);
%      X_init{1} = m{1}*V;
%      [U,V] = pca(m{2},3);
%      X_init{2} = m{2}*V;
%      X_init = [X_init{1} X_init{2}];
    
else
    fprintf('# Initialising X by performing ppca in concatenated observed (scaled) data...\n');
    X_init = ppcaEmbed([m{1} m{2}], latentDimPerModel*2);
end


%-----------------

latentDim = size(X_init,2);

% Free up some memory
clear('Y')



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
    
    %-------- Add dynamics to the model -----
    if dynUsed
        fprintf(1,'# Adding dynamics to the model...\n');
        optionsDyn{i}.type = 'vargpTime';
        optionsDyn{i}.t =timeStampsTraining;
        optionsDyn{i}.inverseWidth=invWidthMultDyn; % Default: 100
        optionsDyn{i}.initX = X_init; % initX; % That's for the dynamics
        
        kern = kernCreate(optionsDyn{i}.t, dynamicKern); % Default: {'rbf','white','bias'}
        
        
        %---- Default values for dynamics kernel
        if strcmp(kern.comp{2}.type, 'white')
            kern.comp{2}.variance = 1e-2; % Usual values: 1e-1, 1e-3
        end
        
        if strcmp(kern.comp{2}.type, 'whitefixed')
            if ~exist('whiteVar')
                whiteVar = 1e-6;
            end
            kern.comp{2}.variance = whiteVar;
            fprintf(1,'# fixedwhite variance: %d\n',whiteVar);
        end
        
        if strcmp(kern.comp{1}.type, 'rbfperiodic')
            if exist('periodicPeriod')
                kern.comp{1}.period = periodicPeriod;
            end
            fprintf(1,'# periodic period: %d\n',kern.comp{1}.period);
        end
        
        % The following is related to the expected number of
        % zero-crossings.(larger inv.width numerator, rougher func)
        if ~strcmp(kern.comp{1}.type,'ou')
            kern.comp{1}.inverseWidth = optionsDyn{i}.inverseWidth./(((max(t{i})-min(t{i}))).^2);
            kern.comp{1}.variance = 1;
        end
        %------
        
        optionsDyn{i}.kern = kern;
        optionsDyn{i}.vardistCovars = vardistCovarsMult; % 0.23 gives true vardist.covars around 0.5 (DEFAULT: 0.23) for the ocean dataset
        
        % Fill in with default values whatever is not already set
        optionsDyn{i} = vargplvmOptionsDyn(optionsDyn{i});
        model{i} = vargplvmAddDynamics(model{i}, 'vargpTime', optionsDyn{i}, optionsDyn{i}.t, 0, 0,optionsDyn{i}.seq);
        
        fprintf(1,'# Further calibration of the initial values...\n');
        model{i} = vargplvmInitDynamics(model{i},optionsDyn{i});
        
        
        model = model{i};
        return 
        %  to also not learn the last kernel's variance
        if numel(kern.comp) > 1 && exist('learnSecondVariance') && ~learnSecondVariance
            fprintf(1,'# The variance for %s in the dynamics is not learned!\n',kern.comp{end}.type)
            model{i}.dynamics.learnSecondVariance = 0;
            model{i}.dynamics.kern.comp{end}.inverseWidth = model{i}.dynamics.kern.comp{1}.inverseWidth/10; %%% TEMP
        end
    end
    
    model{i}.beta=1/(0.01*var(model{i}.m(:)));
    prunedModelInit{i} = vargplvmPruneModel(model{i});
    %disp(model{i}.vardist.covars)
end



%modelInit = model;%%%TEMP

%--  Unify models into a structure
model = svargplvmModelCreate(model);
model.dataSetNames = dataSetNames;
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
%{
fprintf('# Parallel computations w.r.t the submodels!\n');
model.parallel = 1;
model = svargplvmPropagateField(model,'parallel', 1);
%
fprintf('# Parallel computations w.r.t the datapoints!\n');
model.vardist.parallel = 1;
for i=1:model.numModels
    model.comp{i}.vardist.parallel = 1;
end
%}

% Force kernel computations
params = svargplvmExtractParam(model);
model = svargplvmExpandParam(model, params);

%%


display = 1;
%%%% Optimisation
% do not learn beta and sigma_f for few iterations for intitialization
if initVardistIters ~=0
    model.initVardist = 1; model.learnSigmaf = 0;
    model = svargplvmPropagateField(model,'initVardist', model.initVardist);
    model = svargplvmPropagateField(model,'learnSigmaf', model.learnSigmaf);
    fprintf(1,'# Intitiliazing the variational distribution...\n');
    model = svargplvmOptimise(model, display, initVardistIters); % Default: 20
    %fprintf(1,'1/b = %.4d\n',1/model.beta);
    
    modelInitVardist = model;
    model.initVardistIters=initVardistIters;
end

model.initVardist = 0; model.learnSigmaf=1;
model = svargplvmPropagateField(model,'initVardist', model.initVardist);
model = svargplvmPropagateField(model,'learnSigmaf', model.learnSigmaf);


% Optimise the model.
model.iters = 0;
for i=1:length(itNo)
    iters = itNo(i); % default: 2000
    fprintf(1,'\n# Optimising the model for %d iterations (session %d)...\n',iters,i);
    model = svargplvmOptimise(model, display, iters);
    model.iters = model.iters + iters;
    % Save model
    prunedModel = svargplvmPruneModel(model);
    fprintf(1,'# Saving %s\n',fileToSave);
    save(fileToSave, 'prunedModel', 'prunedModelInit');
end

%%

if ~exist('resultsDynamic')  resultsDynamic = 0; end

v = 1;
modelVis = model.comp{v};
if resultsDynamic
    %bar(model.comp{v}.kern.comp{1}.inputScales);
    figure
    % The following causes OUTOFMEMORY exception except from when we prune the
    % video dimensions: (Note: also, lvmVisualise was changed a bit so that no
    % dynamic slides is presented, because otherwise a strange error occurs).
    modelVis.y = Ytr{v};
    [modelP, newHeight, newWidth] = vargplvmReduceVidModel(modelVis, height, width, 2,2);
    lvmVisualise(modelP, [], 'imageVisualise', 'imageModify', [newHeight newWidth],0,0,1);
    clear modelP
end


