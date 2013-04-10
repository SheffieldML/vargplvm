% DEMFINANCE2 Demonstrate Variational GPLVM with dynamics on financial data for multivariate time series (i.e multiple univariate time series)
% DESC The dimensionality of the data is the number of the correlated timeseries. The data number is the points in a timeseries
% (the number of which is the same for all timeseries)
% COPYRIGHT : Andreas C. Damianou, 2011
% VARGPLVM




if ~exist('experimentNo')    experimentNo = 404;    end
if ~exist('itNo')            itNo = 100;            end     % Default: 2000
if ~exist('indPoints')       indPoints = 180;        end     % Default: 180
if ~exist('latentDim')       latentDim = 8;         end     % Default: 8
% Set to 1 to use dynamics or to 0 to use the sta  ndard var-GPLVM
if ~exist('dynUsed')         dynUsed = 1;           end
if ~exist('fixedBetaIters')  fixedBetaIters = 400;    end     % Default: 200
if ~exist('fixInd')          fixInd = 0;            end
if ~exist('dynamicKern')     dynamicKern = {'rbf', 'white', 'bias'}; end
if ~exist('trainModel')      trainModel=1;          end
if ~exist('vardistCovarsMult') vardistCovarsMult=1.5;     end
if ~exist('dataSetName')     dataSetName='USecon';  end
if ~exist('Nstar')           Nstar = 60;            end % Number of test points
if ~exist('makeStationary')  makeStationary = 1;    end
if ~exist('logOfData')       logOfData = 0;         end

% switch dataSetName
%     case 'USecon'
%         load DATA_USecon
%     case 'NYSE2Small'
%         load ../../datasets/finance/nyse_new2Small
% end
Y = lvmLoadData(dataSetName);

%--- FOR DEBUG
if exist('dataToKeep')
    Y = Y(1:dataToKeep,:);
end
%---

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
fprintf(1,'# Iterations (with/without fixed Beta): %d / %s\n',fixedBetaIters,num2str(itNo));
fprintf(1,'# Tie Inducing points: %d\n',fixInd);
fprintf(1,'# Dynamics used: %d\n', dynUsed);
fprintf(1,'# Dataset size used (train/test) : %d / %d \n', size(Ytr,1), size(Yts,1));
if dynUsed
    fprintf(1,'# Dynamics kern: ');
    disp(dynamicKern);
end
fprintf(1,'# Linear correction: %d\n',makeStationary);
fprintf(1,'# Logarithm of tr. data: %d\n',logOfData);
fprintf(1,'#----------------------------------------------------\n');


if (logOfData)
    YtrOriginal = Ytr;
    Ytr = log(Ytr);
end


%----- Remove non-stationarity ----------%
if makeStationary
    t=timeStampsTraining;
    t2=timeStampsTest;
    linearCorrectionTr = zeros(size(Ytr));
    linearCorrectionTs = zeros(size(Yts));
    for q=1:size(Ytr,2)
        % Calculate fit parameters
        [p, ErrorEst] = polyfit(t,Ytr(:,q),2); % Y=p(1)*t^2+p(2)*t+p(3)
        % Evaluate the fit
        curve_fit = polyval(p,t,ErrorEst);
        % Predict (reconstruct) the tr. data itself
        linearCorrectionTr(:,q)=polyval(p,t);
        % Predict for the future (this is needed to "adjust" the GPLVM
        % predictions later)
        linearCorrectionTs(:,q)=polyval(p,t2);
    end
    if ~logOfData
        YtrOriginal = Ytr;
    end
    Ytr = Ytr - linearCorrectionTr;
end
%--------------



%--- Training ---%

options = vargplvmOptions('dtcvar');
options.kern = {'rbfard2', 'bias', 'white'};
options.numActive = indPoints;

%options.scale2var1=1; %%%


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
        optionsDyn.inverseWidth=100; % Usual values: 1-300
        %kern = kernCreate(t, {'rbf', 'white', 'bias'});
        kern=kernCreate(t,dynamicKern);
        kern.comp{2}.variance = 1e-1; % Usual values: 1e-1, 1e-3
        % The following is related to the expected number of zero-crossings.
        kern.comp{1}.inverseWidth = optionsDyn.inverseWidth./(((max(t)-min(t))).^2);
        kern.comp{1}.variance = 1;
        optionsDyn.kern = kern;

        %%%TEMP
        %load TEMP_X_u
        %optionsDyn.X_u = X_u;
        %%%

        optionsDyn.vardistCovars = vardistCovarsMult; % 0.016 gives true vardist.covars around 0.5 (DEFAULT: 0.016), 1.5 gives small covars to avoid noise minimum
        % Set to 0.1 when using MATERN

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
    %modelWriteResult(model, dataSetName, experimentNo);
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
    %     for q=1:model.q, plot(modelInit.vardist.means(:,q)); hold on; plot(model.vardist.means(:,q),'r');
    %         pause(1);
    %         hold off;
    %     end

else % Load an already trained model
    if dynUsed
        %load stickMissing400iters1
        load demUSeconVargplvm404
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

    % REconstruction on tr. data
    VarmuTr = vargplvmPosteriorMeanVar(modelTr, modelTr.vardist.means, modelTr.vardist.covars);
    if makeStationary
        Varmu = Varmu + linearCorrectionTs;
        VarmuTr = VarmuTr + linearCorrectionTr;
    end
    
    if logOfData
        Varmu = exp(Varmu);
        VarmuTr = exp(VarmuTr);
    end

    errsumOnlyTimes = sum((Varmu - YtsOriginal).^2);
    errorOnlyTimes = mean(errsumOnlyTimes);

    % Sampling from a gaussian (this is not correct beause it only accounts for the variance)
    sampMu = Varmu + randn(size(Varmu)).*sqrt(Varsigma);

    % Find the mean sq. error for each datapoint, to see how it increases as
    % we give points further to the future. This time sum accros dimensions
    fprintf(1,'*** Mean errors: %4.6f\n', errorOnlyTimes);
    figure;plot(Varmu-YtsOriginal);
    xlabel('time test datapoints','fontsize',18);
    ylabel('Mean error','fontsize',18);


    %%
    %----------- GARCH comparison -------------------------%
    fprintf(1,'\n# GARCH...\n');
    if makeStationary || logOfData
        Ytr = YtrOriginal;
    end
    predGarch = zeros(size(YtsOriginal));
    predSim = zeros(size(YtsOriginal));
    errorsGarch = zeros(size(YtsOriginal,2),1);
    errorsSim = zeros(size(YtsOriginal,2),1);

    for d=1:size(Y,2)
        YtrRet = price2ret(Ytr(:,d));
        YtsRet = price2ret(YtsOriginal(:,d));

        %Create a specification structure for an ARMA(1,1)/GJR(1,1) model with
        %conditionally t-distributed residuals:
        spec = garchset('VarianceModel','GJR', 'R',1,'M',1,'P',1,'Q',1);
        spec = garchset(spec,'Display','off','Distribution','T');

        %Estimate the parameters of the mean and conditional variance models via
        %garchfit. Make sure that the example returns the estimated residuals and
        %conditional standard deviations inferred from the optimization process,
        %so that they can be used as presample data:
        [coeff,errors,LLF,eFit,sFit] = garchfit(spec,YtrRet);

        horizon = Nstar;  % Define the forecast horizon
        %Call the forecasting engine, garchpred, with the estimated model parameters, coeff, the YtrRet returns, and the forecast horizon:
        [sigmaForecast,meanForecast,sigmaTotal,meanRMSE] = garchpred(coeff,YtrRet,horizon);

        %Simulate 20000 paths (columns):
        nPaths = 20000;  % Define the number of realizations.
        strm = RandStream('mt19937ar','Seed',5489);
        RandStream.setDefaultStream(strm);
        [eSim,sSim,ySim] = garchsim(coeff,horizon,nPaths,[],[],[], eFit,sFit,YtrRet);
        %         pG = ret2price(meanForecast,YtsOriginal(1,d));
        %         pS = ret2price(mean(ySim,2), YtsOriginal(1,d));
        pG = ret2price(meanForecast,Ytr(end,d));
        pS = ret2price(mean(ySim,2), Ytr(end,d));

        predGarch(:,d) = pG(2:end,:);
        predSim(:,d) = pS(2:end,:);

        errorsGarch(d) = mean(sum((predGarch(:,d) - YtsOriginal(:,d)).^2));
        errorsSim(d) = mean(sum((predSim(:,d) - YtsOriginal(:,d)).^2));
        %{
    figure
    plot(t,Ytr(:,d)), hold on
    plot(t2,YtsOriginal(:,d),'c')
    plot(t2,predGarch(:,d),'r')
    plot(t2,predSim(:,d),'g')
    hold off
    errorsGarch(d) = mean(sum((predGarch(:,d) - YtsOriginal(:,d)).^2));
    errorsSim(d) = mean(sum((predSim(:,d) - YtsOriginal(:,d)).^2));
    legend('Ytr','Yts','predGarch','predSim','Location','NorthWest')
    ttl = ['d=' num2str(d) ', errorGarch=', num2str(errorsGarch(d)) ', errorSim=' num2str(errorsSim(d))];
    title(ttl)
    fprintf(1,'# GARCH: %2.4f   |   Sim: %2.4f \n',errorsGarch(d),errorsSim(d));
        %}
    end


    %%

    %---------- Compare with simple linear regression ------------
    fprintf(1,'\n# Linear Regression...\n');
    linRegPlots = 0;

    % If the training data were modified by GPLVM restore them.
    if makeStationary || logOfData
        Ytr = YtrOriginal;
    end
    t=timeStampsTraining;
    t2=timeStampsTest;
    yLinReg = zeros(size(YtsOriginal));
    errLinReg=zeros(size(Ytr,2),1);
    errGPLVM=zeros(size(Ytr,2),1);
    fprintf(1,'# -----------  Error per dimension  --------------\n');
    fprintf(1,'#_dim_|__LR___|_GARCHpred__|__GARCHsim_|__GPLVM___|\n');
    for q=1:size(Ytr,2)

        % Calculate fit parameters
        [p, ErrorEst] = polyfit(t,Ytr(:,q),2); % Y=p(1)*t^2+p(2)*t+p(3)

        % Evaluate the fit
        curve_fit = polyval(p,t,ErrorEst);

        if linRegPlots
            % Plot the data and the fit
            figure
            plot(t,curve_fit,'-',t,Ytr(:,q),'-');

            legend('Polynomial Model','Data','Location','NorthWest');
            xlabel('Time');
            ylabel('TimeSeries Data');
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
        errGPLVM(q)=mean(sum((Varmu(:,q) - YtsOriginal(:,q)).^2));
        fprintf(1,'* d=%d: %4.5f | %4.5f | %4.5f | %4.5f\n', q,errLinReg(q),errorsGarch(q),errorsSim(q), errGPLVM(q));

        if linRegPlots
            figure;plot(yLinReg(:,q)-YtsOriginal(:,q));
            xlabel('time test datapoints','fontsize',18);
            ylabel('Mean error','fontsize',18);
        end


        %%%% SUMMARY
        % Big figure
        scrsz = get(0,'ScreenSize');
        figure('Position',[scrsz(3)/4.86 scrsz(4)/6.666 scrsz(3)/1.6457 scrsz(4)/1.4682])
        if makeStationary
            subplot(2,1,1)
        end
        plot(t,Ytr(:,q)); hold on
        plot(t,VarmuTr(:,q),'r');
        plot(t2,YtsOriginal(:,q),'m');
        plot(t2, Varmu(:,q),'r');
        plot([t;t2], [curve_fit ; yLinReg(:,q) ],'g');
        plot(t2,predGarch(:,q),'k')
        %plot(t2, sampMu(:,q),'c'); %%%%
        legend('Ytr','GPLVMfit','Yts','GPLVM(pred)','LinReg(fit-pred)','GARCH','Location','NorthWest');
        plotTitle=['d = ',num2str(q)];
        title(plotTitle)
        hold off
        if makeStationary
            subplot(2,1,2)
            plot(t,VarmuTr(:,q)-linearCorrectionTr(:,q));
            hold on
            plot(t2,Varmu(:,q)-linearCorrectionTs(:,q),'r');
            title('GPLVM in the stationary space')
        end
        hold off
    end
    fprintf(1,'***MEANS:*** LR: %4.6f | GARCHpred: %4.6f | GARCHsim: %4.6f| GPLVM:%4.6f\n',mean(errLinReg), mean(errorsGarch),mean(errorsSim), mean(errGPLVM));
end


%-----------------




