experiments = [1:5 9:12];
experiments = 13:19
for experimentNo = experiments
    try
        onlyTest = 1;
      %  load(['demOilSvargplvm' num2str(experimentNo) '.mat'])
        dataSetName='oil';
        demClassification
    catch e
        warning(['Could not load experiment ' num2str(experimentNo)])
        getReport(e)
        continue
    end
    svargplvmShowScales(model)
    title(['Experiment: ' num2str(experimentNo)])
    fprintf('# SNR1 = %6f\n', (1/model.comp{1}.beta)/var(model.comp{1}.m(:)));
    fprintf('# SNR2 = %6f\n', (1/model.comp{2}.beta)/var(model.comp{2}.m(:)));
    if isfield(model, 'dynamics') 
        kernDisplay(model.dynamics.kern)
        pause
        figure
        t = model.dynamics.t; Kt = model.dynamics.Kt;X1 = gsamp(zeros(size(Kt,1),1)', Kt,2)';
        plot(X1(:,1), X1(:,2), 'x')
    end
        pause
    close all
    keep('experimentNo', 'experiments')
end

%%

%----------------- Run with itNo


experiments = [1:5 9:12];
for experimentNo = experiments
    try
        onlyTest = 1;
      %  load(['demOilSvargplvm' num2str(experimentNo) '.mat'])
        dataSetName='oil';
        printPlots = 0;
        diaryFile = (['LOGS/LOG_demOilSvargplvm' num2str(experimentNo) '.txt']);
        demClassification
    catch e
        warning(['Could not load experiment ' num2str(experimentNo)])
        getReport(e)
        diary off
        continue
    end
    try
         model = svargplvmReOptimiseModel(model, 0, [600 600]);
    catch e
         warning(['Something went wrong with experiment ' num2str(experimentNo)])
        getReport(e)
        diary off
        continue
    end
    diary off
    keep('experiments','experimentNo')
end


%----------------- Predictions


experiments = [1:5 9:12];
for experimentNo = experiments
    try
        onlyTest = 1;
        doPredictions = 1;
      %  load(['demOilSvargplvm' num2str(experimentNo) '.mat'])
        dataSetName='oil';
        printPlots = 0;
        diaryFile = (['LOGS/LOG_demOilSvargplvm' num2str(experimentNo) '.txt']);
        demClassification
    catch e
        warning(['Could not load experiment ' num2str(experimentNo)])
        getReport(e)
        diary off
        continue
    end
 
    diary off
    keep('experiments','experimentNo')
end

