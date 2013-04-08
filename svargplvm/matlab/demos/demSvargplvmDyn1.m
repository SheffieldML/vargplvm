% DEMSVARGPLVMDYN1 A simple demo on the shared vargplvm with dynamics,
% mostly to check gradients etc.
%
% SHEFFIELDML

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

%------------ Constants -------------
if ~exist('itNo')         ,  itNo = 1;     end     % Default: 2000
if ~exist('indPoints')    ,  indPoints = 5;          end     % Default: 49
if ~exist('initVardistIters'), initVardistIters = 0;      end
if ~exist('mappingKern')   ,  mappingKern = {'rbfard2', 'bias','white'}; end
if ~exist('latentDimPerModel'), latentDimPerModel = 5; end
if ~exist('experimentNo'), experimentNo = 404; end
dynUsed = 1;

%-------------------------------------

dataSetNames = {'Nyse2SmallPart1','Nyse2SmallPart2'};
dataType = 'dyn';
Y = vargplvmLoadData('NYSE2Small');
Y = Y(1:40,:);
half = round(size(Y,1)/2);
Yall{1} = Y(1:half,:);
Yall{2} = Y(half+1:end,:);
t = linspace(0, 2*pi, size(Yall{1}, 1)+1)'; t = t(1:end-1, 1);
timeStampsTraining = t;

clear Y;


numberOfDatasets = length(Yall);

% for i=2:numberOfDatasets
%     if N{i} ~= N{i-1}
%         error('The number of observations in each dataset must be the same!');
%     end
% end

%numberOfDatasets = 1; %%%%%%%%%%%

%%

%-- Options for the models
for i=1:numberOfDatasets
    % Split training and test set
    Ytr{i} = Yall{i}; % no test set
    N{i} = size(Ytr{i},1);
    d{i} = size(Ytr{i},2);
    
    % Set up models
    options{i} = vargplvmOptions('dtcvar');
    options{i}.kern = mappingKern; %{'rbfard2', 'bias', 'white'};
    options{i}.numActive = indPoints;
    options{i}.optimiser = 'scg2';
    options{i}.scaleVal = sqrt(var(Ytr{i}(:)));
    options{i}.enableDgtN = false;
    
    if dynUsed
        optionsDyn{i}.type = 'vargpTime';
        optionsDyn{i}.t=timeStampsTraining;
        optionsDyn{i}.inverseWidth=30;
        %optionsDyn{i}.seq = size(Ytr{i},1);%%%%%%%
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
end

%fprintf('# Initialising X by performing ppca in each observed (scaled) dataset separately and then concatenating...\n');
%X_init{1} = ppcaEmbed(m{1},latentDimPerModel);
%X_init{2} = ppcaEmbed(m{2},latentDimPerModel);
%X_init = [X_init{1} X_init{2}];
X_init = [];
for i=1:numberOfDatasets
    X_init = [X_init ppcaEmbed(m{i}, latentDimPerModel)];
end


%-----------------

latentDim = size(X_init,2);

% Free up some memory
clear('Y')



%-- Create the sub-models.
for i=1:numberOfDatasets
    fprintf(1,'# Creating the model...\n');
    options{i}.initX = X_init;
    model{i} = vargplvmCreate(latentDim, d{i}, Ytr{i}, options{i});
    model{i}.X = X_init;
    model{i} = vargplvmParamInit(model{i}, model{i}.m, model{i}.X);
    model{i}.X = X_init;
    
    inpScales = 5./(((max(model{i}.X)-min(model{i}.X))).^2); % Default 5
    model{i}.kern.comp{1}.inputScales = inpScales;
    
    if strcmp(model{i}.kern.type, 'rbfardjit')
        model{i}.kern.inputScales = model{i}.kern.comp{1}.inputScales;
    end
    params = vargplvmExtractParam(model{i});
    model{i} = vargplvmExpandParam(model{i}, params);
    model{i}.vardist.covars = 0.5*ones(size(model{i}.vardist.covars)) + 0.001*randn(size(model{i}.vardist.covars));
     
    
    model{i}.beta=1/(0.01*var(model{i}.m(:)));
    
    
    
    %-------- Add dynamics to the model -----
    if dynUsed
        model{i} = vargplvmAddDynamics(model{i}, 'vargpTime', optionsDyn{i}, optionsDyn{i}.t, 0, 0,optionsDyn{i}.seq);       
        fprintf(1,'# Further calibration of the initial parameters...\n');
        model{i} = vargplvmInitDynamics(model{i},optionsDyn{i});
    end
    
    
    prunedModelInit{i} = vargplvmPruneModel(model{i});
    %disp(model{i}.vardist.covars)
end



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


% Force kernel computations
params = svargplvmExtractParam(model);
model = svargplvmExpandParam(model, params);

%%

%-------------  Optimisation ---------------
display = 1;
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
    model = svargplvmOptimise(model, display, iters, 'gradcheck', true);
    model.iters = model.iters + iters;
    % Save model
    prunedModel = svargplvmPruneModel(model);
    fprintf(1,'# Saving %s\n',fileToSave);
    save(fileToSave, 'prunedModel', 'prunedModelInit');
end

