
allDataToKeep = [100:100:1000];
experimentNo = 5000;
for jj = 1:length(allDataToKeep)
    randn('seed', 1e4);
    rand('seed', 1e4);
    labelErrorsAll = 0;
    labelErrorsNNAll = 0;
    cmat_relAll = zeros(3,3);
    dataToKeep = allDataToKeep(jj);
    disp('###########-------------------------------------------------------############');
    disp(['################  dataToKeep: ' num2str(dataToKeep) ' ###################']);
    try
 	diary(['demOilSvargplvm' num2str(experimentNo) '.txt']);
        testDataPerClass = -1;
        totalTrials = 1;
	    displayIters = true;
        indPoints = 120;
        dataSetName = 'oil';
        baseKern = {'rbfardjit','rbfardjit'};
        initX = 'pca';
        initial_X = 'concatenated';
        initSNR = 100;
	if dataToKeep < 30
		latentDim = 5;
	else
		latentDim = 6;
	end
        initVardistIters = 75 + min(dataToKeep,80);
        itNo = 250 + min(dataToKeep,500);
        scale2var1 = 1;
        doPredictions = true;
        % 'labels' can be given as the constrain type. In that case, during
        % test time an uninformative prior must be used.
        dynamicsConstrainType = {};
        demClassificationSvargplvm3

        cmat_relAll = cmat_relAll ./ repmat(sum(cmat_relAll, 2), 1, 3);
        disp('*** Conf. Matrix all:')
        cmat_relAll
        fprintf(' *** meanLabelErrors:   %.1f\n *** meanLabelErrorsNN: %.1f\n\n', labelErrorsAll/totalTrials, labelErrorsNNAll/totalTrials)

	experimentNo = experimentNo + 1;
    catch e
        warning(['!! Problem with experiment' num2str(experimentNo)])
        getReport(e)
    end
    diary off
    keep('jj','allDataToKeep', 'experimentNo');
end



%%
subsets = allDataToKeep;

 
svargplvmErrors = [];
NNerrors = [];
varianceAll = [];

allExperiments = [5000:5009];
for i=1:length(allExperiments)
    try
        curExp = allExperiments(i);
        load(['demOilSvargplvmPred' num2str(curExp)])
        
        % fprintf('Exp: %d. ErrorGplvm:%.3f - ErrorGplvmLogs: %.3f\n',curExp, errors.gplvm, svargplvmErrors(i));
        svargplvmErrors = [svargplvmErrors errors.gplvm];
        NNerrors = [NNerrors errors.NN];
        curVar = sum(sum(predictions.varX))/length(predictions.varX(:));
        varianceAll = [varianceAll curVar];
    end
end


accuracySvargplvm = (1000 - svargplvmErrors)/ length(allExperiments);
accuracyNN = (1000 - NNerrors)/ length(allExperiments);

plot(subsets, accuracySvargplvm, '.--', 'MarkerSize', 18)
hold on
plot(subsets, accuracyNN, '.--r','MarkerSize', 18)
hold off
set(gca,'FontSize',12)
legend('MRD','NN','FontSize',12)
xlabel('Number of training points','FontSize',12)
ylabel('Accuracy %','FontSize',12)
