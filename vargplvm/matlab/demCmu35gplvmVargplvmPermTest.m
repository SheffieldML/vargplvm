% DEMCMU35GPLVMVARGPLVMPERMTEST Run variational GPLVM with dynamics on CMU35 data - permutation test to test the bound.
%
% COPYRIGHT :  Andreas C. Damianou, Michalis K. Titsias, 2011
%
% SEEALSO : demCmu35vargplvmReconstructTaylor.m
% VARGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

% Define constants (in a manner that allows other scripts to parametrize
% this one).
if ~exist('experimentNo') ,  experimentNo = 404;      end
if ~exist('itNo')         ,  itNo = [200 1000 1300];              end     % Default: 2000
if ~exist('indPoints')    ,  indPoints = 100;          end     % Default:
if ~exist('latentDim')    ,  latentDim = 9;          end
% Set to 1 to use dynamics or to 0 to use the standard var-GPLVM
if ~exist('dynUsed')      ,  dynUsed = 1;             end

if ~exist('fixedBetaIters'), fixedBetaIters = 50;      end     % DEFAULT:
% Set to 1 to tie the inducing points with the latent vars. X
if ~exist('fixInd')        ,     fixInd = 0;                             end
if ~exist('baseKern')   ,     baseKern = {'rbfard2', 'white'}; end
if ~exist('dynamicKern')   ,     dynamicKern = {'rbf', 'bias', 'white'}; end
if ~exist('reconstrIters') ,     reconstrIters = 2000;                   end
if ~exist('learnVariance'),     learnVariance =0;   end
%if ~exist('vardistCovarsMult'),  vardistCovarsMult=1.3;                  end
if ~exist('initX'),     initX ='ppca';   end
if ~exist('doReconstr'),     doReconstr=1;   end
if ~exist('permTestMult'), permTestMult = 5; end



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
%fprintf(1,'# VardistCovarsMult: %d\n',vardistCovarsMult);
fprintf(1,'# Reconstruction iterations: %d\n', reconstrIters);
fprintf(1,'# Dataset size used (train) : %d \n', size(Y,1));
fprintf(1,'# Tie Inducing points: %d\n',fixInd);
fprintf(1,'# InitX: %s\n',initX);
fprintf(1,'# Base kern: ');
disp(baseKern);
fprintf(1,'# Dynamics used: %d\n', dynUsed);
if dynUsed
    fprintf(1,'# Dynamics kern: ');
    disp(dynamicKern);
end

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
timeStampsTraining = []; timeStampsTrainingOrig = [];
dt=0.05;
for i=1:length(seq)
    tOrig = ([0:(seq(i)-prevSeq)].*dt)';
    timeStampsTrainingOrig = [timeStampsTrainingOrig ;tOrig];
    
    t = 0;
   % Add gaussian noise to the timestep so that time only goes forward but
   % not necessarily in equal steps (check diff(timeStampsTraining))
    for j=2:seq(i)-prevSeq + 1
        t = [t; t(end) + dt + abs(randn)*dt*permTestMult];
    end
    timeStampsTraining = [timeStampsTraining ;t];
    
    prevSeq = seq(i)+1;

end
timeStampsTest = ([0:size(Ytest,1)-1].*dt)';


% Set up model
options = vargplvmOptions('dtcvar');
%options.kern = {'rbfard2', 'bias', 'white'};
options.kern = baseKern;
options.numActive = indPoints; % DEfault: 100

options.optimiser = 'scg';

d = size(Y, 2);
fprintf(1,'# Creating the model...\n');
if fixInd
    options.fixInducing=1;
    options.fixIndices=1:indPoints;
end
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
    %   optionsDyn.vardistCovars = vardistCovarsMult;
    optionsDyn.seq = seq;
    optionsDyn.learnVariance = learnVariance;
    optionsDyn.initX = initX;
    optionsDyn.regularizeMeans = 0;
    
    % Dynamic kernel:
    kern = kernCreate(t, dynamicKern); % Default: {'rbf','white','bias'}
    
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
        kern.comp{1}.inverseWidth = optionsDyn.inverseWidth./(((max(t)-min(t))).^2);
        kern.comp{1}.variance = 1;
    end
    optionsDyn.kern = kern;
    
    if exist('vardistCovarsMult')
        optionsDyn.vardistCovars = vardistCovarsMult;
    end
    
    
    
    % Fill in with default values whatever is not already set
    optionsDyn = vargplvmOptionsDyn(optionsDyn);
    model = vargplvmAddDynamics(model, 'vargpTime', optionsDyn, optionsDyn.t, 0, 0,optionsDyn.seq);
    
    fprintf(1,'# Further calibration of the initial parameters...\n');
    model = vargplvmInitDynamics(model,optionsDyn);
end
modelInit = model;


% do not learn beta for few iterations for intitilization
if fixedBetaIters > 0
    model.learnBeta = 0;
    display = 1;
    fprintf(1,'# Intitiliazing the model (fixed beta) %d iterations...\n',fixedBetaIters);
    model = vargplvmOptimise(model, display, fixedBetaIters);
    model.fixedBetaIters = fixedBetaIters;
    disp('# Saving model after optimising beta...')
    try
        vargplvmWriteResult(model, model.type, dataSetName, experimentNo);
    catch e
        warning('Error while saving model: ')
        fprintf('%s\n',e.message)
    end
end

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
    try
        vargplvmWriteResult(model, model.type, dataSetName, experimentNo);
    catch e
        warning('Error while saving model: ')
        fprintf('%s\n',e.message)
    end
end

%bar(model.kern.comp{1}.inputScales)
%prefix = 'scales';
%saveAllOpenFigures(['Results/CMU/NEW/' num2str(experimentNo) '/'], prefix,1)


% Reconstruction can also be done separately calling demCmu35vargplvmReconstructTaylor
if doReconstr
    %---- RECONSTRUCTION ---%
    save 'TEMPExperimentNo.mat' 'experimentNo'
    clear
    load 'TEMPExperimentNo'
    delete TEMPExperimentNo.mat
    fprintf('# Taylor Reconstruction for expNo:%d\n', experimentNo);
    demCmu35vargplvmReconstructTaylor
    %prefix = [num2str(experimentNo) '/'];
    %saveAllOpenFigures('Results/CMU/NEW/', prefix,1)
end
