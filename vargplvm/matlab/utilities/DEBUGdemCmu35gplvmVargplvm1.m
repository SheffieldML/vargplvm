% DEMSTICKVARGPLVM1 Run variational GPLVM on CMU35 data.

% SHEFFIELDML

clear
% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dynUsed=1;
iters = 50; % Default: 1000

% Get the sequence numbers.
[Y, lbls] = lvmLoadData('cmu35WalkJog');
seq = cumsum(sum(lbls)) - [1:31];

dataSetName = 'cmu35gplvm';
experimentNo = 100;

% load data
[Y, lbls, Ytest, lblstest] = lvmLoadData(dataSetName);

% Scaling the data so that they have zero mean and unit var.
%disp('# Scaling data...')
% The following 2 lines are done by taylor
% origBias = mean(Y);
% origScale = 1./sqrt(var(Y));
%origBias = mean(Y);
%origScale = ones(1, size(Y,2));

% !!!! DO NOT CHANGE THE ORIGINAL DATA !!!!
% If YOU WANT TO USE TAYLOR SCALING YOU CAN DO EITHER AFTER modelCreate, BY DEFINING THERE THE NEW 
% model.scale and model.bias AND RE-DEFINING model.m FROM THE ORIGINAL model.y (***which  is never changed***))
% (THIS WILL USE TAYLOR'S NORMALIZATION FOR TRAINING) 
% OR
% AT REDICTION TIME BY TRANFORMING THE FINAL PREDICTION OF THE MODEL AND
% THE CORRESPONDING TEST DATA
% (THIS ALLOWS YOU TO USE ANY NORMALIZATION YOU PREFER FOR TRAINING) 

%scale = ones(size(scale));
%Y = Y - repmat(origBias, size(Y, 1), 1);
%Ytest = Ytest - repmat(origBias, size(Ytest, 1), 1);
%Y = Y.*repmat(origScale, size(Y, 1), 1);
%Ytest = Ytest.*repmat(origScale, size(Ytest, 1), 1);
%Y=Y(1:190,:); %%%%% temp: for gradchek
%seq=seq(1:2);%%%% temp: for gradchek


%%%% Remove some sequences
seqFrom=1;
seqEnd=3;


if seqFrom ~= 1
    Yfrom = seq(seqFrom-1)+1;
else
    Yfrom = 1;
end
Yend=seq(seqEnd);
Y=Y(Yfrom:Yend,:);
seq=seq(seqFrom:seqEnd);
seq=seq - ones(1,length(seq)).*(Yfrom-1);


% Fix times:
prevSeq = 1;
timeStampsTraining = [];
dt=0.05;
for i=1:length(seq)
    t = ([0:(seq(i)-prevSeq)].*dt)';
    prevSeq = seq(i)+1;
    timeStampsTraining = [timeStampsTraining ;t];
end;
timeStampsTest = ([0:size(Ytest,1)-1].*dt)';

% Set up model
options = vargplvmOptions('dtcvar');
options.kern = {'rbfard2', 'white'};
options.numActive = 30; % DEfault: 100
% model.m (the actual data for the model) is transformed to have zero mean and unit var.
options.scale2var1 = 1;

options.optimiser = 'scg';
latentDim = 5; % Try for more, e.g. 10
d = size(Y, 2);
fprintf(1,'# Creating the model...\n');
model = vargplvmCreate(latentDim, d, Y, options);
model = vargplvmParamInit(model, model.m, model.X);
model.beta=100;
model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));
if dynUsed
    model = addDefaultVargpTimeDynamics(model, timeStampsTraining,seq);
    %model = addDefaultVargpTimeDynamics(model, timeStampsTraining);%%%temp: for gradchek
end


if dynUsed
    disp('# Further calibration of the parameters...')
    model.dynamics.kern.comp{1}.inverseWidth = 30./(((max(model.dynamics.t)-min(model.dynamics.t))).^2);
    params = vargplvmExtractParam(model);
    model = vargplvmExpandParam(model, params);
    % Initialize barmu
    initFunc = str2func([options.initX 'Embed']);
    X = initFunc(model.m, model.q);
    vX = var(X);
    for q=1:model.q
        Lkt = chol(model.dynamics.Kt + 0.01*vX(q)*eye(model.N))';
        model.dynamics.vardist.means(:,q) = Lkt'\(Lkt\X(:,q));
      %  model.dynamics.vardist.means(:,q) = randn(model.N,1);  %Lkt'\(Lkt\X(:,q));  %%%%% This is probably worse
    end
    % smaller lengthscales
    model.kern.comp{1}.inputScales = 5./(((max(X)-min(X))).^2);
    params = vargplvmExtractParam(model);
    model = vargplvmExpandParam(model, params);
    % inducing point need to initilize based on model.vardist.means
    perm = randperm(model.k);
    model.X_u = model.vardist.means(perm(1:model.k),:);
    params = vargplvmExtractParam(model);
    model = vargplvmExpandParam(model, params);
end
modelInit = model;

fprintf(1,' # LatentDim: %d\n # Ind.Pts: %d\n # ItersNoBeta: %d\n',latentDim, options.numActive);

% Optimise the model.
display = 1;
model.learnBeta = 1;
fprintf(1,'# Optimising the model for %d iterations...\n',iters);
model = vargplvmOptimise(model, display, iters);
fprintf('1 over b=%.4d\n', 1/model.beta);
% Save the results.
fprintf(1,'# Saving model after doing %d iterations',iters)
modelWriteResult(model, dataSetName, experimentNo);
