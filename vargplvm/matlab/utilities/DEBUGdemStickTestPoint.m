% DEMSTICKVARGPLVM1 Run variational GPLVM on stick man data.

% SHEFFIELDML

clear; 

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'stick';

clear timeStamps; % in case it's left from a previous experiment

% load data
[Y, lbls] = lvmLoadData(dataSetName);
dims = size(Y,2);
N=size(Y,1);

% Declare constants (in a manner that allows other scripts to parametrize
% this one.
if ~exist('expNo')
    expNo = 100;
end
if ~exist('itNo')
    itNo = 4; % Default: 2000
end
if ~exist('indPoints')
    indPoints = 50; % Default: 50
end
if ~exist('latentDim')
    latentDim = 10;
end
if ~exist('dynUsed')
    dynUsed =1;
end
if ~exist('dataToKeep')
    dataToKeep = size(Y,1);
end
if dataToKeep == -1
    dataToKeep = size(Y,1);
end

experimentNo = expNo;
printDiagram = 0;
displayDynamic = 0;
trainModel = 1;% Set it to 1 to retrain the model. Set it to 0 to load an already trained one.

% % Temp: for fewer data
% if dataToKeep < size(Y,1)
%     fprintf(1,'# Using only a subset of %d datapoints...\n', dataToKeep);
%     Y = Y(1:dataToKeep,:);
% end

% Set training and test sets
Ytr = Y(1:end-5,:);
Yts = Y(end-4:end,:);


t = linspace(0, 2*pi, size(Y, 1)+1)';  
t = t(1:end-1, 1);

timeStampsTraining = t(1:end-5,1);
timeStampsTest = t(end-4:end,1);


% Set up model
options = vargplvmOptions('dtcvar');
options.kern = {'rbfard2', 'bias', 'white'};
options.numActive = indPoints; 
%options.tieParam = 'tied';  

options.optimiser = 'scg';
%latentDim = 10;
d = size(Y, 2);



 

    % demo using the variational inference method for the gplvm model
    fprintf(1,'# Creating the model...\n');
    model = vargplvmCreate(latentDim, d, Ytr, options);
    %
    model = vargplvmParamInit(model, model.m, model.X); 
    model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));


    % Stickman data are equally-spaced (t(n+1) = t(n) + 0.0083).
    %-------- Add dynamics to the model -----
    if dynUsed
        model = addDefaultVargpTimeDynamics(model, timeStampsTraining);
    end
    %------------------



    if ~isempty(model.dynamics)
        % A better way to initialize the  kernel hyperparameter,
        % especially lengthscale, should be to learn it by running few iterations of
        % GP regression marginal likelihood maximization given as data the PCA output
        % (the default value is jsut reasonable and it sets the inverse lenthscale to quite 
        % small value so the GP dynamic prior is weaker (less smoothed) at
        % initializations
        model.dynamics.kern.comp{1}.inverseWidth = 200./(((max(model.dynamics.t)-min(model.dynamics.t))).^2);
        params = vargplvmExtractParam(model);
        model = vargplvmExpandParam(model, params);    
        % Initialize barmu     
        initFunc = str2func([options.initX 'Embed']);
        X = initFunc(model.m, model.q);
        vX = var(X); 
        for q=1:model.q
            Lkt = chol(model.dynamics.Kt + 0.01*vX(q)*eye(model.N))';    
            % barmu = inv(Kt + s*I)*X, so that  mu = Kt*barmu =  Kt*inv(Kt +
            % s*I)*X, which is jsut the GP prediction, is temporally smoothed version 
            % of the PCA latent variables X (for the data Y)
            model.dynamics.vardist.means(:,q) = Lkt'\(Lkt\X(:,q));    
        end
        % smaller lengthscales
        model.kern.comp{1}.inputScales = 5./(((max(X)-min(X))).^2);
        params = vargplvmExtractParam(model);
        model = vargplvmExpandParam(model, params);
        % inducing point need to initilize based on model.vardist.means
        perm = randperm(model.k); 
        model.X_u = model.vardist.means(perm(1:model.k),:);
    
        % Note from Michalis: It seems that optimization gets stuck to the trivial local minima where
        % everything is expained by the noise in the likelihood (ie. 1/model.beta 
        % becomes equal to var(model.Y(:))). With fixed model.beta, 
        % you can escape from this local minima, but it not clear how to 
        % automatically initialize so that to resolve this problem 
        params = vargplvmExtractParam(model);
        model = vargplvmExpandParam(model, params);
    end
    modelInit = model;
    
    
%     % do not learn beta for few iterations for intitilization
%     model.learnBeta = 0;
%     display = 1;
%     fprintf(1,'# Intitiliazing the model (fixed beta)...\n');
%     model = vargplvmOptimise(model, display, 200);
    
    % Optimise the model.
    iters = itNo; % default: 2000
    display = 1;
    model.learnBeta = 1;
    fprintf(1,'# Optimising the model...\n');
    model = vargplvmOptimise(model, display, iters);
    
    % fprintf(1,'Plotting inferred latent GPs (intial with blue and final  with red)...\n');
    % for q=1:model.q, plot(modelInit.vardist.means(:,q)); hold on; plot(model.vardist.means(:,q),'r'); 
    %     pause(1); 
    %     hold off; 
    % end
    



% Save the results.
fprintf(1,'# Saving the model...\n');
%modelWriteResult(model, dataSetName, experimentNo);
 capName = dataSetName;
 capName(1) = upper(capName(1));
 modelType = model.type;
 modelType(1) = upper(modelType(1));
 save(['dem' capName modelType num2str(experimentNo) '.mat'], 'model');


if exist('printDiagram') & printDiagram
  lvmPrintPlot(model, lbls, dataSetName, experimentNo);
end

if exist('displayDynamic') & displayDynamic
  % load connectivity matrix
  [void, connect] = mocapLoadTextData('run1');
  % Load the results and display dynamically.
  lvmResultsDynamic(model.type, dataSetName, experimentNo, 'stick', connect)
end

