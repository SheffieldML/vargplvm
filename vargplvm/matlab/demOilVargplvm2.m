% DEMOILVARGPLVM2 Run variational GPLVM on oil data.

% VARGPLVM

% Fix seeds
randn('seed', 1e6);
rand('seed', 1e6);

dataSetName = 'oil';
experimentNo = 2;
printDiagram = 1;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

%%% TEMP: for fewer data
Y = Y(1:100,:); 

% training and test sets
Ntr = 0.7*floor(size(Y,1)); % 70% of the data for training  
perm = randperm(size(Y,1)); 
Ytr = Y(perm(1:Ntr),:);      lblsTr = lbls(perm(1:Ntr),:);
Yts = Y(perm(Ntr+1:end),:);  lblsTs = lbls(perm(Ntr+1:end),:);


% Set up model
options = vargplvmOptions('dtcvar');
options.kern = {'rbfard2', 'bias', 'white'};
options.numActive = 20; % Default: 50

options.optimiser = 'scg';
latentDim = 10;
d = size(Y, 2);

model = vargplvmCreate(latentDim, d, Ytr, options);
%
model = vargplvmParamInit(model, model.m, model.X); 
model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));

iters = 10; % Default: 2000
display = 1;

model = vargplvmOptimise(model, display, iters);

iters = 100;
display = 0;

% 50% missing outputs from the each test point
numIndPresent = round(0.5*model.d); 
% patrial reconstruction of test points
for i=1:size(Yts,1)
    %
    % randomly choose which outputs are present
    permi = randperm(model.d);
    indexPresent =  permi(1:numIndPresent);
    indexP(i,:) = indexPresent;
    
    % initialize the latent point using the nearest neighbour 
    % from he training data
    dst = dist2(Yts(i,indexPresent), Ytr(:,indexPresent));
    [mind, mini] = min(dst);
    
    Init(i,:) = model.vardist.means(mini,:);
    % create the variational distribtion for the test latent point
    vardistx = vardistCreate(model.vardist.means(mini,:), model.q, 'gaussian');
    vardistx.covars = 0.2*ones(size(vardistx.covars));
   
    % optimize mean and vars of the latent point 
    model.vardistx = vardistx;
    vardistx = vargplvmOptimisePoint(model, vardistx, Yts(i,indexPresent), indexPresent, display, iters);
    Testmeans(i,:) = vardistx.means; 
    Testcovars(i,:) = vardistx.covars;
    
    % reconstruct the missing outputs  
    [mu, sigma] = vargplvmPosteriorMeanVar(model, vardistx);
    Varmu(i,:) = mu; 
    Varsigma(i,:) = sigma; 
    %
end


capName = dataSetName;;
capName(1) = upper(capName(1));
modelType = model.type;
modelType(1) = upper(modelType(1));
save(['dem' capName modelType num2str(experimentNo) '.mat'], 'model', 'perm', 'indexP', 'Varmu', 'Varsigma');

    
