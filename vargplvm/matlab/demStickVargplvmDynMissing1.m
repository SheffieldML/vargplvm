% DEMSTICKVARGPLVMDYNMISSING1 Run variational GPLVM on stick man data with
% the option to add dynamics. Predict on test data from which some
% dimensions are partly observed. Evaluate the reconstruction.

% VARGPLVM


%clear;

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'stick';

clear timeStamps; % in case it's left from a previous experiment

% load data
[Y, lbls] = lvmLoadData(dataSetName);
dims = size(Y,2);

% Define constants (in a manner that allows other scripts to parametrize
% this one).
if ~exist('experimentNo')    experimentNo = 4044;    end
if ~exist('itNo')            itNo = 50;            end     % Default: 2000
if ~exist('indPoints')       indPoints = 10;        end     % Default: 50
if ~exist('latentDim')       latentDim = 8;        end     % Default: 10
% Set to 1 to use dynamics or to 0 to use the sta  ndard var-GPLVM
if ~exist('dynUsed')         dynUsed = 1;           end
if ~exist('fixedBetaIters')  fixedBetaIters = 50;   end     % Default: 200
if ~exist('reconstrIters')   reconstrIters = 50;    end      %Default: 200
if ~exist('fixInd')          fixInd = 0;            end    


%%%%%%%%%%% TEMP
%delete TEMPInducingMeansDist.mat
%delete TEMPInducingMeansDistNorm.mat
%delete TEMPInducingIndices.mat
%delete TEMPLikelihoodTrace.mat
%%%%%%%%%%%%%%%%%%%

printDiagram = 0;
displayDynamic = 0;
trainModel = 1;  % Set it to 1 to retrain the model. Set it to 0 to load an already trained one.


% Set training and test sets
Nstar = 5; % number of test points
overL = 0; % allow overlap (0 means no overlap)
Ytr = Y(1:end-Nstar,:);
Yts = Y(end-(Nstar-1+overL):end,:);

% Corresponding timestamps (artificial and equally spaced for this demo)
% The real timestamps have t(n+1) = t(n) + 0.0083
t = linspace(0, 2*pi, size(Y, 1)+1)';
t = t(1:end-1, 1);
timeStampsTraining = t(1:end-Nstar,1);
timeStampsTest = t(end-(Nstar-1+overL):end,1);


fprintf(1,'\n#----------------------------------------------------\n');
fprintf(1,'# Dataset: %s\n',dataSetName);
fprintf(1,'# ExperimentNo: %d\n', experimentNo);
fprintf(1,'# Inducing points: %d\n',indPoints);
fprintf(1,'# Latent dimensions: %d\n',latentDim);
fprintf(1,'# Iterations (with/without fixed Beta): %d / %d\n',fixedBetaIters,itNo);
fprintf(1,'# Tie Inducing points: %d\n',fixInd);
fprintf(1,'# Dynamics used: %d\n', dynUsed);
fprintf(1,'# Dataset size used (train/test) : %d / %d \n', size(Ytr,1), size(Yts,1));
fprintf(1,'#----------------------------------------------------\n');


% Set up model
options = vargplvmOptions('dtcvar');
options.kern = {'rbfard2', 'bias', 'white'};
options.numActive = indPoints;

options.optimiser = 'scg';
d = size(Y, 2);

if trainModel
    fprintf(1,'# Creating the model...\n');
    if fixInd
        options.fixInducing=1;
        options.fixIndices=1:size(Ytr,1);
    end
    model = vargplvmCreate(latentDim, d, Ytr, options);
    model = vargplvmParamInit(model, model.m, model.X);
    model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));
    
    %-------- Add dynamics to the model -----
    if dynUsed
        optionsDyn.type = 'vargpTime';
        optionsDyn.t=timeStampsTraining;
        optionsDyn.inverseWidth=200;
        
        % Fill in with default values whatever is not already set
        optionsDyn = vargplvmOptionsDyn(optionsDyn);
        model = vargplvmAddDynamics(model, 'vargpTime', optionsDyn, optionsDyn.t, 0, 0,optionsDyn.seq);
        
        fprintf(1,'# Further calibration of the initial parameters...\n');
        model = vargplvmInitDynamics(model,optionsDyn);
    end
    model.beta=1/(0.01*var(model.m(:)));
    modelInit = model;
    
    % It seems that optimization gets stuck to the trivial local minima where
    % everything is expained by the noise in the likelihood (ie. 1/model.beta
    % becomes equal to var(model.Y(:))). With fixed model.beta,
    % we can escape from this local minima, but it not clear how to
    % automatically initialize so that to resolve this problem
    
    % do not learn beta for few iterations for intitilization
    model.learnBeta = 0;
    display = 1;
	if (fixedBetaIters ~= 0)
        fprintf(1,'# Intitiliazing the model (fixed beta)...\n');
        model = vargplvmOptimise(model, display, fixedBetaIters);
    end
    
    % Optimise the model.
    model.learnBeta = 1;
    iters = itNo;
    fprintf(1,'# Optimising the model...\n');
    model = vargplvmOptimise(model, display, iters);
    
    % fprintf(1,'Plotting inferred latent GPs (intial with blue and final  with red)...\n');
    % for q=1:model.q, plot(modelInit.vardist.means(:,q)); hold on; plot(model.vardist.means(:,q),'r');
    %     pause(1);
    %     hold off;
    % end
    
else % Load an already trained model
    if dynUsed
        %load stickMissing400iters1
        %load  stickMissing1iter
        load demStickVargplvm404
    else
        load demStickVargplvm101
    end
end

if ~dynUsed %%% TEMP!!
    if trainModel
        % Save the results.
        fprintf(1,'# Saving the model...\n');
        %modelWriteResult(model, dataSetName, experimentNo);
        capName = dataSetName;
        capName(1) = upper(capName(1));
        modelType = model.type;
        modelType(1) = upper(modelType(1));
        save(['dem' capName modelType num2str(experimentNo) '.mat'], 'model');
    end
    if exist('printDiagram') & printDiagram
        lvmPrintPlot(model, lbls, dataSetName, experimentNo);
    end
    displayDynamic=1; printDiagram=1; %%%
    if exist('displayDynamic') & displayDynamic
        % load connectivity matrix
        [void, connect] = mocapLoadTextData('run1');
        % Load the results and display dynamically.
        lvmResultsDynamic(model.type, dataSetName, experimentNo, 'stick', connect)
    end
    return
end


return %%%%


%%%%%%%%%%%%%------ Reconstruction --------%%%%%%%%%%

display = 1;
iters=reconstrIters;
fprintf(1, '# Removing some outputs randomly from test data...\n');


% 50% missing outputs from each test point
numIndPresent = round(0.5*model.d);
indexP = [];
Init = [];
YtsOriginal = Yts;

fprintf(1, '# Partial reconstruction of test points...\n');
% randomly choose which outputs are present
permi = randperm(model.d);
indexPresent =  permi(1:numIndPresent);
%indexP(i,:) = indexPresent;

for i=1:size(Yts,1)
    % initialize the latent points using the nearest neighbour
    % from he training data
    dst = dist2(Yts(i,indexPresent), Ytr(:,indexPresent));
    [mind, mini(i)] = min(dst);
end

% create the variational distribtion for the test latent point
if dynUsed
    vardistx = vardistCreate(model.dynamics.vardist.means(mini,:), model.q, 'gaussian');
    vardistx.covars = model.dynamics.vardist.covars(mini,:);%0.2*ones(size(vardistx.covars));
else
    % Not tested yet!!!!
    vardistx = vardistCreate(model.vardist.means(mini,:), model.q, 'gaussian');
    vardistx.covars = model.vardist.covars(mini,:);%0.2*ones(size(vardistx.covars));
end
% optimize mean and vars of the latent point
model.vardistx = vardistx;
indexMissing = setdiff(1:model.d, indexPresent);

Yts(:,indexMissing) = NaN;
if dynUsed
    model.dynamics.t_star = timeStampsTest;
    % The returned quantities are not the original S and mu, but the
    % new parameters mu_bar and lambda.
    modelTr = model;
    [x, varx, model] = vargplvmOptimisePoint(model, vardistx, Yts, display, iters);
    
    % keep the optmized parameters
    barmu = x;
    lambda = varx;
    
    % Get the variational means and variances for the new test data and
    % update the model to be prepared for prediction
    [x, varx, modelUpdated] = vargplvmDynamicsUpdateModelTestVar(model, barmu, lambda);
    
else % static vargplvm
    [x, varx] = vargplvmOptimisePoint(model, vardistx, Yts, display, iters);
end


%%%%%%%%-------- Evaluate reconstruction ------------%%%%%%%%%
% Find the square error
Testmeans = x;
Testcovars = varx;
[mu, sigma] = vargplvmPosteriorMeanVar(modelUpdated, x, varx);
Varmu = mu;
Varsigma = sigma;
errsumFull = sum((Varmu(:,indexMissing) - YtsOriginal(:,indexMissing)).^2);
errorFull = mean(errsumFull);

% Compute also the error without updating the model after the re-training with
% partially observed test data (recall that this affects the variational
% in the training data)
[mu, sigma] = vargplvmPosteriorMeanVar(modelTr, x, varx);
Varmu1 = mu;
Varsigma1 = sigma;
errsumNoUpdate = sum((Varmu1(:,indexMissing) - YtsOriginal(:,indexMissing)).^2);
errorNoUpadate = mean(errsumNoUpdate);

% Baseline prediction using only the information in the test time points
[Testmeans2 Testcovars2] = vargplvmPredictPoint(modelTr.dynamics, modelTr.dynamics.t_star);
[Varmu2, Varsigma2] = vargplvmPosteriorMeanVar(modelTr, Testmeans2, Testcovars2);
errsumOnlyTimes = sum((Varmu2(:,indexMissing) - YtsOriginal(:,indexMissing)).^2);
errorOnlyTimes = mean(errsumOnlyTimes);


% Find the mean sq. error for each datapoint, to see how it increases as
% we give points further to the future. This time sum accros dimensions
errSq = sum((Varmu(:,indexMissing) - YtsOriginal(:,indexMissing)).^2 ,2);
errPt = errSq ./ size(Varmu,2);
fprintf(1,'*** Mean errors: Full=%4.6f, No updating after testing=%4.6f, Baseline=%4.6f\n', errorFull, errorNoUpadate, errorOnlyTimes);
figure;plot(errPt);
xlabel('time test datapoints','fontsize',18);
ylabel('Mean error','fontsize',18);


% find the largest dimensions
[max, imax] = sort(model.kern.comp{1}.inputScales,'descend');

N = size(model.X,1);
Nstar = size(Yts,1);
figure;
fprintf(1,'The two larger latent dims after re-training with partial test data. Red are the test points\n');
plot(modelUpdated.X(1:model.N,imax(1)), modelUpdated.X(1:model.N,imax(2)),'--rs', 'Color','b');
hold on
plot(Testmeans(:,imax(1)),Testmeans(:,imax(2)),'--rs', 'Color','r');
title('Visualization of the latent space after re-training with partial test data');
hold off

figure;
fprintf(1,'The two larger latent dims. Green are baseline test predictions.\n');
plot(model.X(1:model.N,imax(1)), model.X(1:model.N,imax(2)),'--rs', 'Color','b')
hold on
plot(Testmeans2(:,imax(1)),Testmeans2(:,imax(2)),'--rs', 'Color','g');
title('Visualization of the latent space without re-training');
hold off;

fprintf(1,'# Reconstructing the data and showing error bars with dashed lines...\n');
cnt = 0;
figure
for ii=1:size(indexMissing,2)
    i = indexMissing(ii);
    if mod(ii,16)==0, figure;  cnt = 1;
    else, cnt = cnt + 1;
    end
    subplot(4,4,cnt);
    plot(YtsOriginal(:,i), 'Color', 'r'); hold on;
    plot(Varmu(:,i));  plot(Varmu(:,i) - 2*sqrt(Varsigma(:,i)),'b:'); plot(Varmu(:,i) + 2*sqrt(Varsigma(:,i)),'b:');
    plot(Varmu2(:,i),'g');  plot(Varmu2(:,i) - 2*sqrt(Varsigma2(:,i)),'g:'); plot(Varmu2(:,i) + 2*sqrt(Varsigma2(:,i)),'g:');
    titlestring = ['Missing output:' num2str(i)];
    title(titlestring);
    axis tight;
    %legend('Original data','Reconstructed data','Error bars');
    pause(0.5);
    hold off;
end

%%%%%%%%%% end: evaluate reconstruction

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
