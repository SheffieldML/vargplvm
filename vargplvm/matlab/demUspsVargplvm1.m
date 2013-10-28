% DEMUSPSVARGPLVM1 Demonstrate variational GPLVM on USPS data. 

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
