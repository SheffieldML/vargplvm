
clear
randn('seed', 1e5);
rand('seed', 1e5);


experimentNo = 27; % Like 5 but random subset
diaryFile = (['LOGS/LOG_demOilSvargplvm' num2str(experimentNo) '.txt']);
delete(diaryFile)
diary(diaryFile)
clear diaryFile
      
        
totalTrials = 15; % 10
for trialNo = 1:totalTrials
    disp(['############################  TRIAL NO ' num2str(trialNo) ' #####################'])
    disp('')
    try
        dataSetName = 'oil';
        baseKern = {'rbfardjit','rbfardjit'};
        dataToKeep = 50;
        initX = 'pca';
        initial_X = 'concatenated';
        initSNR = 100;
        latentDim = 6;
        initVardistIters = 530; % 500
        itNo = 1000; %1000
        testDataPerClass = 100; % 100
        scale2var1 = 1;
        doPredictions = true;
        dynamicsConstrainType = {}; %%% dynUsed = 0
        tic; demClassification2; toc
        close all
    catch e
        warning(['!! Problem with experiment' num2str(experimentNo)])
        getReport(e)
    end
    
    scales1All{trialNo} = model.comp{1}.kern.inputScales;
    scales2All{trialNo} = model.comp{2}.kern.inputScales;
    
    
    indsAll{trialNo} = indTrial;
    labelErrorsAll{trialNo} = labelErrors;
    labelSharingErrorsAll{trialNo} = labelSharingErrors;
    labelErrorsNNAll{trialNo} = labelErrorsNN;
    wrongIndicesAll{trialNo} = wrongIndices;
    wrongIndicesNNAll{trialNo} = wrongIndicesNN;
    SNR1{trialNo} = (1/model.comp{1}.beta)/var(model.comp{1}.m(:));
    SNR2{trialNo} = (1/model.comp{2}.beta)/var(model.comp{2}.m(:));
    cmat_relAll{trialNo} = cmat_rel;
    keep('experimentNo','totalTrials','trialNo', 'indsAll', 'labelErrorsAll', ...
        'labelSharingErrorsAll', 'labelErrorsNNAll', 'SNR1', 'SNR2', ...
        'scales1All', 'scales2All', 'wrongIndicesAll', 'wrongIndicesNNAll', 'cmat_relAll');
end

%{
for i=1:totalTrials
    totalPoints = length(indsAll{i});
    class1 = length(find(labelsTest(indsAll{i}) == 1)) * (100 / totalPoints);
    class2 = length(find(labelsTest(indsAll{i}) == 2)) * (100 / totalPoints);
    class3 = length(find(labelsTest(indsAll{i}) == 3)) * (100 / totalPoints);
    fprintf('# TrialNo = %d: This dataset has class 1,2,3 in percentage: (%.1f, %.1f, %.1f) percent. \n', i, class1, class2, class3);
end
%}

save(['scriptOilRandom' num2str(experimentNo) '.mat']);
diary off