% DEMCMU35GPLVMVARGPLVM1 Run variational GPLVM with dynamics on CMU35 data.
%
% COPYRIGHT :  Andreas C. Damianou, Michalis K. Titsias, 2011
%
% SEEALSO : demCmu35gplvmVargplvm3.m
% VARGPLVM






clear
% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

%---- Experiment parameters -----
dynUsed=1;
indPoints=100;
latentDim=9;
experimentNo = 1;
fixedBetaIters = 50;
reconstrIters = 2000;
itNo = [200 1000 1300];
%----------------------------------

% Get the sequence numbers.
[Y, lbls] = lvmLoadData('cmu35WalkJog');
seq = cumsum(sum(lbls)) - [1:31];

dataSetName = 'cmu35gplvm';


% load data
[Y, lbls, Ytest, lblstest] = lvmLoadData(dataSetName);


fprintf(1,'\n#----------------------------------------------------\n');
fprintf(1,'# Dataset: %s\n',dataSetName);
fprintf(1,'# ExperimentNo: %d\n', experimentNo);
fprintf(1,'# Inducing points: %d\n',indPoints);
fprintf(1,'# Latent dimensions: %d\n',latentDim);
fprintf(1,'# Iterations (with/without fixed Beta): %d / %s\n',fixedBetaIters, num2str(itNo));
fprintf(1,'# Reconstruction iterations: %d\n', reconstrIters);
fprintf(1,'# Dataset size used (train) : %d \n', size(Y,1));
fprintf(1,'#----------------------------------------------------\n');





% %%% Uncomment this code to get a subset of the sequences
% seqFrom=2;
% seqEnd=4;
%
% if seqFrom ~= 1
%     Yfrom = seq(seqFrom-1)+1;
% else
%     Yfrom = 1;
% end
% Yend=seq(seqEnd);
% Y=Y(Yfrom:Yend,:);
% seq=seq(seqFrom:seqEnd);
% seq=seq - ones(1,length(seq)).*(Yfrom-1);



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
%options.kern = {'rbfard2', 'bias', 'white'};
options.kern = {'rbfard2', 'white'};
options.numActive = indPoints; % DEfault: 100

options.optimiser = 'scg';

d = size(Y, 2);
fprintf(1,'# Creating the model...\n');
model = vargplvmCreate(latentDim, d, Y, options);
model = vargplvmParamInit(model, model.m, model.X);
model.beta=1/(0.01*var(model.m(:)));
model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));

model.reconstrIters = reconstrIters;

%-------- Add dynamics to the model -----
if dynUsed
    optionsDyn.type = 'vargpTime';
    optionsDyn.t=timeStampsTraining;
    optionsDyn.inverseWidth=30;
    optionsDyn.seq = seq;
    
    % Fill in with default values whatever is not already set
    optionsDyn = vargplvmOptionsDyn(optionsDyn);
    model = vargplvmAddDynamics(model, 'vargpTime', optionsDyn, optionsDyn.t, 0, 0,optionsDyn.seq);
    
    fprintf(1,'# Further calibration of the initial parameters...\n');
    model = vargplvmInitDynamics(model,optionsDyn);
end
modelInit = model;


% do not learn beta for few iterations for intitilization
model.learnBeta = 0;
display = 1;
fprintf(1,'# Intitiliazing the model (fixed beta) %d iterations...\n',fixedBetaIters);
model = vargplvmOptimise(model, display, fixedBetaIters);
model.fixedBetaIters = fixedBetaIters;
disp('# Saving model after optimising beta...')
modelWriteResult(model, dataSetName, experimentNo);

% Optimise the model.
display = 1;
model.learnBeta = 1;
model.iters=0;

for i=1:length(itNo)
    iters = itNo(i); % Default: 1000
    fprintf(1,'# Optimising the model for %d iterations (session %d)...\n',iters,i);
    model = vargplvmOptimise(model, display, iters);
    fprintf('1/b=%.4d var(model.m)=%.4d\n', 1/model.beta, var(model.m(:)));
    model.iters = model.iters + iters;
    model.date = date;
    % Save the results.
    fprintf(1,'# Saving model after doing %d iterations',iters)
    modelWriteResult(model, dataSetName, experimentNo);
end
%{
iters = 150; % Default: 1000
fprintf(1,'# Optimising the model for %d iterations...\n',iters);
model = vargplvmOptimise(model, display, iters);
fprintf('1/b=%.4d var(model.m)=%.4d\n', 1/model.beta, var(model.m(:)));
model.iters = model.iters + iters;
% Save the results.
fprintf(1,'# Saving model after doing %d iterations',iters)
modelWriteResult(model, dataSetName, experimentNo);

iters = 600; % Default: 1000
fprintf(1,'# Optimising the model for %d iterations...\n',iters);
model = vargplvmOptimise(model, display, iters);
fprintf('1/b=%.4d var(model.m)=%.4d\n', 1/model.beta, var(model.m(:)));
model.iters = model.iters + iters;
% Save the results.
fprintf(1,'# Saving model %s after doing %d iterations',dataSetName,iters)
modelWriteResult(model, dataSetName, experimentNo);

iters = 1400; % Default: 1000
fprintf(1,'# Optimising the model for %d iterations...\n',iters);
model = vargplvmOptimise(model, display, iters);
fprintf('1/b=%.4d var(model.m)=%.4d\n', 1/model.beta, var(model.m(:)));
% Save the results.
fprintf(1,'# Saving model %s after doing %d iterations',dataSetName,iters)
modelWriteResult(model, dataSetName, experimentNo);
%}
