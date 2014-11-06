% DEMUSPSVARGPLVM1 Demonstrate variational GPLVM on USPS data. 
% Copyright: Michalis Titsias, Neil Lawrence, Andreas Damianou, 2010 - 2014
% VARGPLVM

% Fix seeds
rng(1e5, 'v4');

dataSetName = 'usps';
experimentNo = 1;
printDiagram = 1;

% load data
[YTrain, lblsTrain, YTest, lblsTest] = lvmLoadData(dataSetName);

nClasses = size(lblsTrain,2);
% Set up model
options = vargplvmOptions('dtcvar');
options.kern = {'rbfard2', 'white'};
options.numActive = 50; 
%options.tieParam = 'tied';

options.optimiser = 'scg';
latentDim = 10;
d = size(YTrain, 2);

capName = dataSetName;
capName(1) = upper(capName(1));
modelType = 'vargplvm'; % varmodel{1}.type;
modelType(1) = upper(modelType(1));
modelFile = ['dem' capName modelType num2str(experimentNo) '.mat'];

if exist(modelFile, 'file') % Training data exists, just load it.
    load(modelFile);
else % Do the training:
    iters = 1000;
    display = 1;
    
    varmodel = cell(nClasses,1);
    % create a separate vargplvm for each digit
    for i=1:nClasses
        %
        fprintf('\n\n# Training for Class # %d\n\n',i)
        Y = YTrain(lblsTrain(:,i)==1,:);
        
        model = vargplvmCreate(latentDim, d, Y, options);
        %
        model = vargplvmParamInit(model, model.m, model.X);
        
        model = vargplvmOptimise(model, display, iters);
        
        varmodel{i} = model;
        %
    end
end

iters = 100;
display = 0;

save(modelFile, 'varmodel');

% measure performance on test data 
prob = zeros(size(YTest,1),nClasses);
TestError = 0;

% New and faster way: give the whole test matrix at once, rather than using
% a for loop through each test element. This is almost 10 times faster then
% the old way (using 2 Matlab workers on a 4-core machine).

% Compute the approximate class conditional density for each digit
tic
for i=1:nClasses
    sprintf('--------\nComputing test probability for class=%d\n', i);
    prob(:,i) = vargplvmProbabilityCompute(varmodel{i}, YTest, 0, iters);
end
[maxP C] = max(prob,[],2);
for n=1:size(YTest,1)
    TestError = TestError + ~lblsTest(n,C(n));
end
fprintf('TestError=%d,\tTime=%fseconds.\n', TestError, toc);
%
save(['dem' capName modelType num2str(experimentNo) '.mat'], 'varmodel', 'prob', 'TestError');


%% Nearest Neighbour baseline
dists = dist2(YTrain,YTest);
[~,positions] = min(dists);
Ypred_NN = lblsTrain(positions,:);
NNError = 0;
for n=1:size(YTest,1)
    pp = find(Ypred_NN(n,:));
    NNError = NNError + ~lblsTest(n, pp);
end

%% Training of logistic regression classifier
labelsTrain = transformLabels(lblsTrain)';
labelsTest = transformLabels(lblsTest)';
Nstar = size(YTest,1);

for i=1:nClasses
    fprintf('\n # LogReg training for class # %d\n', i)
    lb = zeros(size(YTrain,1),1);
    lb(labelsTrain == i) = 1;
    B{i} = glmfit(YTrain, lb,'binomial','logit'); % Logistic regression
end

% Prediction of each binary classifier
Ypred_logReg = zeros(size(lblsTest));
for i=1:nClasses
    Ypred_logReg(:,i) = glmval(B{i},YTest,'logit')';
end

% Replace predictions with maximum probability (ie, make a decision)
[~,ind]=max(Ypred_logReg');
LogRegError = 0;
for i=1:size(YTest,1)
    LogRegError = LogRegError + (ind(i) ~= labelsTest(i));
end

% Print all results
fprintf('Var-GPLVM = %d (%f %%) \n', TestError, TestError*100/Nstar);
fprintf('NN        = %d mistakes (%f %% misclassified) \n', NNError, NNError*100/Nstar);
fprintf('Log. Reg. = %d mistakes (%f %% misclassified) \n', LogRegError, LogRegError*100/Nstar);
