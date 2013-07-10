% DEMROBOTWIRELESSVARGPLVMDYN1 Run variational GPLVM on robot wireless
% data.
%
% DESC: Run variational GPLVM on robot wireless data with the option to
% also add dynamics to the model.
%
% COPYRIGHT : Michalis K. Titsias, 2009-2011
% COPYRIGHT : Neil D. Lawrence, 2009-2011
% COPYRIGHT : Andreas C. Damianou, 2010-2011
%
% VARGPLVM


%clear;


% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'robotWireless';
clear timeStamps; % in case it's left from a previous experiment

%-- Constants
if ~exist('trainModel'), trainModel=1; end
if ~exist('indPoints'), indPoints = 70; end % Default: 50
if ~exist('latentDim'), latentDim = 10; end
if ~exist('dynUsed'),  dynUsed = 1; end
if ~exist('printDiagram'), printDiagram = 1; end
if ~exist('experimentNo'), experimentNo = 404; end
if ~exist('doPredictions'), doPredictions = 0; end
if ~exist('itNo'), itNo = 2000; end
if ~exist('initVardistIters'), initVardistIters = 180; end
if ~exist('dynamicKern')   ,     dynamicKern = {'rbf', 'white', 'bias'}; end
if ~exist('vardistCovarsMult'), vardistCovarsMult = 2.5; end

% load data
%[Y, lbls] = lvmLoadData(dataSetName);
%%% Like lvmLoadData but also parse times.
fprintf(1,'# Preparing the dataset...\n');
baseDir = datasetsDirectory;
dirSep = filesep;
[Ydat, timeStampsdat, wireless_x, wireless_y, storedMacs] = parseWirelessData([baseDir 'uw-floor.txt']);
Ydat = (Ydat + 85)/15;
lbls = 'connect';

if ~exist('dataToKeep')
    dataToKeep = size(Ydat,1); % Default should be: 215
end

Y = Ydat(1:dataToKeep, :);
timeStampsTraining = timeStampsdat(1:dataToKeep);
Ytest = Ydat((dataToKeep+1):end, :);
timeStampsTest = timeStampsdat((dataToKeep+1):end);
%%%


% Set up model
options = vargplvmOptions('dtcvar');
options.kern = {'rbfard2', 'bias', 'white'};
options.numActive = indPoints; % Default: 50

options.optimiser = 'scg';
d = size(Y, 2);

capName = dataSetName;
capName(1) = upper(capName(1));
modelType = 'Vargplvm';
saveName = ['dem' capName modelType num2str(experimentNo) '.mat'];


if trainModel
    % demo using the variational inference method for the gplvm model
    fprintf(1,'# Creating the model...\n');
    model = vargplvmCreate(latentDim, d, Y, options);
    %
    model = vargplvmParamInit(model, model.m, model.X);
    model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));
    if dynUsed
        optionsDyn.type = 'vargpTime';
        optionsDyn.t=timeStampsTraining;
        optionsDyn.inverseWidth=20;
        optionsDyn.kern = dynamicKern;
        optionsDyn.vardistCovars = vardistCovarsMult;
        
        % Fill in with default values whatever is not already set
        optionsDyn = vargplvmOptionsDyn(optionsDyn);
        model = vargplvmAddDynamics(model, 'vargpTime', optionsDyn, optionsDyn.t, 0, 0,optionsDyn.seq);
        
        fprintf(1,'# Further calibration of the initial parameters...\n');
        model = vargplvmInitDynamics(model,optionsDyn);
    end
    modelInit = model;
    
    
    %     for q=1:model.q, plot(X(:,q)); hold on; plot(model.vardist.means(:,q),'r');
    %         plot(model.X_u(:,q),'+r');  pause(1); hold off;
    %     end
    %
    
    %model.dynamics.vardist.covars = ones(size(model.dynamics.vardist.covars)); % Good initialization
    
    % Optimise the model.
    model.initVardist = 1; model.learnSigmaf = 0;
    display = 1;
    fprintf(1,'# Optimising the model (initialising var.distr) for %d iters...\n',initVardistIters);
    model = vargplvmOptimise(model, display, initVardistIters);
    
    model.initVardist = 0; model.learnSigmaf = 1;
    
    fprintf(1,'# Optimising the model for %d iters...\n',itNo);
    model = vargplvmOptimise(model, display, itNo);
    
    fprintf(1,'# Saving the model...(%s)\n',saveName);
    
    prunedModel = vargplvmPruneModel(model);
    save(saveName, 'prunedModel');
    
    
    % order wrt to the inputScales
    mm = vargplvmReduceModel(model,2);
    % plot the two largest twe latent dimensions
    if exist('printDiagram') & printDiagram
        lvmPrintPlot(mm, lbls, capName, experimentNo);
    end
    
else
    load(saveName);
    model = vargplvmRestorePrunedModel(prunedModel, Y);
end



% fprintf(1,'Plotting inferred latent GPs (intial with blue and final  with red)...\n');
% for q=1:model.q, plot(X(:,q)); hold on; plot(model.vardist.means(:,q),'r');
%     plot(model.X_u(:,q),'+r');  pause(1); hold off;
% end

%%%% Reconstruct training data %%%%%%%%
% mu  = vargplvmPosteriorMeanVar(model, model.vardist.means, model.vardist.covars);
% Varmutr=mu;
% % TODO

if ~(doPredictions && dynUsed)
    return
end

%%%% predictions
Yts = Ytest;
Varmu = [];
Varsigma = [];

fprintf(1, '# Prediction...\n');
t_star = timeStampsTest;
[x varx] = vargplvmPredictPoint(model.dynamics, t_star);


% find the largest dimensions
[max, imax] = sort(model.kern.comp{1}.inputScales,'descend');

N = size(model.X,1);
Nstar = size(Yts,1);
figure;
fprintf(1,'The two larger latent dims after re-training with partial test data. Red are the test points\n');
plot(model.X(1:model.N,imax(1)), model.X(1:model.N,imax(2)),'--rs', 'Color','b');
hold on
plot(x(:,imax(1)),x(:,imax(2)),'--rs', 'Color','r');
title('Visualization of the latent space after re-training with partial test data');
hold off

figure;

% for i=1:size(t_star,1)
%     [mu, sigma] = vargplvmPosteriorMeanVar(model, x(i,:), varx(i,:));
%     Varmu(i,:) = mu;
%     Varsigma(i,:) = sigma;
%     %
% end

[mu, sigma] = vargplvmPosteriorMeanVar(model, x, varx);
Varmu = mu;
Varsigma = sigma;

% Find the square error
errsum = sum((Varmu - Yts).^2);
% Devide by the total number of test points to find a mean
error = errsum / size(Varmu,1);

% Find the mean sq. error for each datapoint, to see how it increases as
% we give points further to the future. This time sum accros dimensions
errSq = sum((Varmu - Yts).^2 ,2);
errPt = errSq ./ size(Varmu,2);

fprintf(1,'*** Mean error: %d\n', mean(error));

% See how the latent space looks like
N = size(model.X,1);
Nstar = size(t_star,1);
newX = zeros( N+Nstar, latentDim );
newX(1:N,:) = model.X(:,:);
newX(N+1:N+Nstar, :) = x(:,:);
figure;plot(errPt);
xlabel('time test datapoints','fontsize',18);
ylabel('Mean error','fontsize',18);
%---


% Now we will plot the latent space INCLUDING the newly predicted points.
model.X = newX;
modelTEMP = vargplvmCreate(latentDim, d, Ydat, options);
model.m = modelTEMP.m;
model.y = modelTEMP.y;
model.N = modelTEMP.N;
model.vardist.means = [model.vardist.means; x];
model.vardist.covars = [model.vardist.covars; varx];
%%%%%%%%% end: predictions



% order wrt to the inputScales
mm = vargplvmReduceModel(model,2);
% plot the two largest twe latent dimensions
if exist('printDiagram') & printDiagram
    lvmPrintPlot(mm, lbls, capName, experimentNo);
    %lvmResultsDynamic(model.type, dataSetName, experimentNo, 'robotWireless','vector')
end

