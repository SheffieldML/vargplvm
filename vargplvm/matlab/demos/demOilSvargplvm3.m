% DEMOILSVARGPLVM3 Classification demo on the oil data for various subsets
% of the data and comparison with NN results.
%
% COPYRIGHT: Andreas C. Damianou, 2013
%
% VARGPLVM

allDataToKeep = [70 100:100:1000];
experimentNo = 6000;
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
        delete(['demOilSvargplvm' num2str(experimentNo) '.txt']);
        diary(['demOilSvargplvm' num2str(experimentNo) '.txt']);
        testDataPerClass = -1;
        totalTrials = 1;
	    displayIters = true;
        indPoints = 120;
        dataSetName = 'oil';
        baseKern = {'rbfardjit','rbfardjit'};
        initX = 'ppca';
        initial_X = 'custom';
        initSNR = 100;
	if dataToKeep < 30
        latentDimPerModel = {3,0};
	else
		latentDimPerModel = {4,0};
	end
        initVardistIters = 75 + min(dataToKeep,80);
        itNo = 300 + min(dataToKeep,500);
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
    %prunedModel = svargplvmPruneModel(model);
    %save[];
    diary off
    keep('jj','allDataToKeep', 'experimentNo');
end



%%
%subsets = allDataToKeep;
subsets = [];
 
svargplvmErrors = [];
NNerrors = [];
varianceAll = [];
scalesAll = {};

allExperiments = [6000:6012];
for i=1:length(allExperiments)
    try
        curExp = allExperiments(i);
        load(['demOilSvargplvmPred' num2str(curExp)])
        
        % fprintf('Exp: %d. ErrorGplvm:%.3f - ErrorGplvmLogs: %.3f\n',curExp, errors.gplvm, svargplvmErrors(i));
        scalesAll{end+1} = scales;
        svargplvmErrors = [svargplvmErrors errors.gplvm];
        NNerrors = [NNerrors errors.NN];
        curVar = sum(sum(predictions.varX))/length(predictions.varX(:));
        varianceAll = [varianceAll curVar];
        subsets = [subsets allDataToKeep(i)];
    catch e
        fprintf('%s\n',e.message)
    end
end


accuracySvargplvm = (1000 - svargplvmErrors) ./ 10;
accuracyNN = (1000 - NNerrors)/ 10;

plot(subsets, accuracySvargplvm, '.--', 'MarkerSize', 18)
hold on
plot(subsets, accuracyNN, '.--r','MarkerSize', 18)
hold off
set(gca,'FontSize',12)
legend('MRD','NN','FontSize',12)
xlabel('Number of training points','FontSize',12)
ylabel('Accuracy %','FontSize',12)

%{
for i=1:length(scalesAll)
    subplot(1,2,1); bar(scalesAll{i}{1}); subplot(1,2,2); bar(scalesAll{i}{2}); pause; 
end
%}