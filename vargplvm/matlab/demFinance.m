% DEMFINANCE Demonstrate variational GPLVM with dynamics on financial data for a univariate time series (one time series)
% DESC This demo is not really of interest the dimensionality of the data is D=1
% COPYRIGHT : Andreas C. Damianou, 2011
% VARGPLVM

load garchdata
dataSetName='garchdata';

% Y=price2ret(NASDAQ); % This is the return
Y=NASDAQ; % This is the price series

clear DEM2GBP
clear NYSE
clear NASDAQ



if ~exist('experimentNo')    experimentNo = 404;    end
if ~exist('itNo')            itNo = 100;            end     % Default: 2000
if ~exist('indPoints')       indPoints = 80;        end     % Default: 50
if ~exist('latentDim')       latentDim = 1;        end     % Default: 10
% Set to 1 to use dynamics or to 0 to use the sta  ndard var-GPLVM
if ~exist('dynUsed')         dynUsed = 1;           end
if ~exist('fixedBetaIters')  fixedBetaIters = 0;   end     % Default: 200
if ~exist('fixInd')          fixInd = 0;            end    
if ~exist('trainModel')      trainModel=1;          end
if ~exist('dynamicKern')     dynamicKern = {'rbf', 'white', 'bias'}; end
if ~exist('dataToKeep')      dataToKeep=500;         end
if ~exist('Nstar')           Nstar=30;         end
if ~exist('kernInvWidth')    kernInvWidth=100; end
%Y=Y(1:dataToKeep);
Y=Y(200:500); %%%%TEMP!!!


Nstar = 30; % number of test points
overL = 0; % allow overlap (0 means no overlap)
Ytr = Y(1:end-Nstar,:);
Yts = Y(end-(Nstar-1+overL):end,:);

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
if dynUsed
    fprintf(1,'# Dynamics kern: ');
    disp(dynamicKern);
    fprintf(1,'# kernInvWidth: %d\n', kernInvWidth);
end
fprintf(1,'# Dataset size used (train/test) : %d / %d \n', size(Ytr,1), size(Yts,1));
fprintf(1,'#----------------------------------------------------\n');


options = vargplvmOptions('dtcvar');
options.kern = {'rbfard2', 'bias', 'white'};
options.numActive = indPoints;

options.optimiser = 'scg';
d = size(Y, 2);

if trainModel
    fprintf(1,'# Creating the model...\n');
    options.scaleVal = sqrt(var(Ytr(:)));

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
        optionsDyn.inverseWidth=kernInvWidth; % The bigger this is, the rougher the function
        
        kern=kernCreate(t,dynamicKern);
        kern.comp{2}.variance = 1e-1; % Usual values: 1e-1, 1e-3
        % The following is related to the expected number of zero-crossings.
        kern.comp{1}.inverseWidth = optionsDyn.inverseWidth./(((max(t)-min(t))).^2);
        kern.comp{1}.variance = 1;
        optionsDyn.kern = kern;

        
        optionsDyn.vardistCovars = 0.2; % For 500 dataToKeep 0.2 gives initial vardist.covars closer to 0 rather to 1, helps escaping noise minimum
        

        
        % Fill in with default values whatever is not already set
        optionsDyn = vargplvmOptionsDyn(optionsDyn);
        model = vargplvmAddDynamics(model, 'vargpTime', optionsDyn, optionsDyn.t, 0, 0,optionsDyn.seq);
        
        fprintf(1,'# Further calibration of the initial parameters...\n');
        model = vargplvmInitDynamics(model,optionsDyn);
    end
    model.beta=1/(0.01*var(model.m(:)));
    modelInit = model;
    
        fprintf(1,'# Model initial vardist. covars:\n');
    disp(model.vardist.covars)
    
   
    %--
    capName = dataSetName;
    capName(1) = upper(capName(1));
    modelType = model.type;
    modelType(1) = upper(modelType(1));
    fileName=['dem' capName modelType num2str(experimentNo) '.mat'];
    %--
    
    % do not learn beta for few iterations for intitilization
    model.learnBeta = 0;
    display = 1;
	if (fixedBetaIters ~= 0)
        fprintf(1,'# Intitiliazing the model (fixed beta) for %d iters...\n', fixedBetaIters);
        model = vargplvmOptimise(model, display, fixedBetaIters);
    end
    model.fixedBetaIters=fixedBetaIters;
    
    model.iters=0;
    
    % Optimise the model.
    model.learnBeta = 1;
    for i=1:length(itNo)
        iters=itNo(i);
        fprintf(1,'# Optimising the model for %d iters...(session %d)\n', iters,i);
        model = vargplvmOptimise(model, display, iters);
        model.iters = model.iters+iters;
        fprintf(1,'# Saving the model (%s)...\n', fileName );
        save(fileName, 'model', 'modelInit');
    end
    
    
    % fprintf(1,'Plotting inferred latent GPs (intial with blue and final  with red)...\n');
    % for q=1:model.q, plot(modelInit.vardist.means(:,q)); hold on; plot(model.vardist.means(:,q),'r');
    %     pause(1);
    %     hold off;
    % end
    
else % Load an already trained model
    if dynUsed
        %load stickMissing400iters1
        %load  stickMissing1iter
    else
        %load stickMissing400iters1
    end
end

fprintf(1,'# Optimisation complete. \n# 1/b=%.5d \n# var(model.m(:))=%.5d\n',1/model.beta, var(model.m(:)))


%----- Test (only for the dynamic case)
if dynUsed
    YtsOriginal = Yts;
    model.dynamics.t_star = timeStampsTest;
    % The returned quantities are not the original S and mu, but the
    % new parameters mu_bar and lambda.
    modelTr = model;
    [Testmeans Testcovars] = vargplvmPredictPoint(modelTr.dynamics, modelTr.dynamics.t_star);
    [Varmu, Varsigma] = vargplvmPosteriorMeanVar(modelTr, Testmeans, Testcovars);
    
    errsumOnlyTimes = sum((Varmu - YtsOriginal).^2);
    errorOnlyTimes = mean(errsumOnlyTimes);
    
    
    % Find the mean sq. error for each datapoint, to see how it increases as
    % we give points further to the future. This time sum accros dimensions
    
    fprintf(1,'*** Mean errors: %4.6f\n', errorOnlyTimes);
    figure;plot(Varmu-YtsOriginal);
    xlabel('time test datapoints','fontsize',18);
    ylabel('Mean error','fontsize',18);
    
    
    %---------- Compare with simple linear regression ------------
    linRegPlots = 1;
    
    t=timeStampsTraining;
    t2=timeStampsTest;
    yLinReg = zeros(size(YtsOriginal));
    errLinReg=zeros(size(Ytr,2),1);
    fprintf(1,'# Error per dimension, GPLVM vs Lin. Reg\n');
    for q=1:size(Ytr,2)
        
        % Calculate fit parameters
        [p, ErrorEst] = polyfit(t,Ytr(:,q),2); % Y=p(1)*t^2+p(2)*t+p(3)
        
        % Evaluate the fit
        curve_fit = polyval(p,t,ErrorEst);
        
        if linRegPlots
            % Plot the data and the fit
            figure
            plot(t,curve_fit,'-',t,Ytr(:,q),'-');
        end
        
        % Predict for the future
        yLinReg(:,q)=polyval(p,t2);
        errsumLinReg = sum((yLinReg(:,q) - YtsOriginal(:,q)).^2);
        errLinReg(q) = mean(errsumLinReg);
        
        if linRegPlots
            hold on;
            plot(t2, yLinReg(:,q),'-');
            plot(t2,YtsOriginal(:,q),'-','Color','r');
            hold off;
            
            legend('Polynomial Model','Data','Location','NorthWest');
            xlabel('Time');
            ylabel('TimeSeries Data');
        end
        
        % Find the mean sq. error for each datapoint, to see how it increases as
        % we give points further to the future. This time sum accros dimensions
        
        fprintf(1,'***Dim %d: ErrorLR: %4.6f   |    ErrorGPLVM=%4.6f\n', q,errLinReg(q),mean(sum((Varmu(:,q) - YtsOriginal(:,q)).^2)));
        
        if linRegPlots
            figure;plot(yLinReg(:,q)-YtsOriginal(:,q));
            xlabel('time test datapoints','fontsize',18);
            ylabel('Mean error','fontsize',18);
        end
    end
end

mu=model.vardist.means;
% Rescale the mean
mu = mu.*repmat(model.scale, size(model.vardist.means,1), 1);
% Add the bias back in
mu = mu + repmat(model.bias, size(model.vardist.means,1), 1);


plot(timeStampsTraining, Ytr,'b'); hold on; 
plot(timeStampsTraining, mu,'r');
plot(timeStampsTest,Varmu, 'r'); 
plot(timeStampsTraining, curve_fit,'m');
plot(timeStampsTest, yLinReg,'m');
plot(timeStampsTest, YtsOriginal, 'g');
legend('Ytr','GPLVM(fit)','GPLVM(pred)','Lin.Reg.(fit)','Lin.Reg(pred)','YtsOriginal','Location','SouthEast');
title('Summary')
hold off

