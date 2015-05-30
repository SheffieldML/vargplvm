% DEMCLASSIFICATIONSVARGPLVM3
% COPYRIGHT Andreas C. Damianou, 2012
% SEEALSO: demClassificationSvargplvm2, demClassificationSvargplvm,  demClassificationSvargplvm4
% VARGPLVM

% Like demClassificationSvargplvm2 but we train
% with RANDOM subsets of the training data what necessarily have the same number
% of instances per class.

if ~exist('doPredictions'),    doPredictions = true; end
if ~exist('testDataPerClass'), testDataPerClass = -1; end %%% All data
if ~exist('printPlots'), printPlots = false; end
if ~exist('totalTrials'), totalTrials = 10; end
if ~exist('displayIters'), displayIters = false; end


%addpath(genpath('../../bc-vargplvm/matlab'));
%addpath(genpath('../../vargplvm/matlab'));

% This script initialises the options structure 'globalOpt'.
svargplvm_init;

globalOpt.dataSetNames = {'oilData', 'oilLabels'};
globalOpt.dataSetName = 'oil';


[YtrFull, lblsTrainFull] = svargplvmLoadData('oil');
labelsTrainFull = transformLabels(lblsTrainFull);
[YtestFull, lblsTest] = svargplvmLoadData('oilTest');
labelsTest = transformLabels(lblsTest);

if isfield(globalOpt, 'normaliseData') && globalOpt.normaliseData
    warning('normalising data independently!')
    YtrFull = utils_normaliseData(YtrFull);
    YtestFull = utils_normaliseData(YtestFull);
end

if globalOpt.dataToKeep == size(YtrFull,1)
    globalOpt.dataToKeep = -1;
end


if testDataPerClass ~= -1
[curIndTest, Yts{1}, labelsTest, Yts{2}] = utils_subsetOfClasses(YtestFull, ...
    testDataPerClass,labelsTest, lblsTest);
else
    Yts{1} = YtestFull;
    Yts{2} = lblsTest;
end
lblsTest = Yts{2};

Yts{2}(Yts{2} == 0) = -1;

for trialNo = 1:totalTrials
    if globalOpt.dataToKeep ~= -1
        curData = {'dataToKeep',globalOpt.dataToKeep};
    else
        curData = globalOpt.dataPerClass;
    end

    disp(['################  TRIAL NO ' num2str(trialNo) ' ###################']);
    if (globalOpt.dataToKeep == -1 && globalOpt.dataPerClass == -1) %%%
        Yall{1} = YtrFull; Yall{2} = lblsTrainFull;  %%%
        labelsTrain = labelsTrainFull; %%%
    else
        [curIndAll{trialNo}, Yall{1}, labelsTrain, Yall{2}] = utils_subsetOfClasses(YtrFull, ...
        curData,labelsTrainFull, lblsTrainFull);
    end
    lblsTrain = Yall{2};
    % Transform e.g. 1 0 0 to 1 -1 -1
    Yall{2}(Yall{2} == 0) = -1;


    %---------------------------------------------------------------


    numberOfDatasets = length(Yall);
    globalOpt.indPoints = min(globalOpt.indPoints, size(Yall{1},1));

    %-- Load datasets
    for i=1:numberOfDatasets
        Y = Yall{i};
        dims{i} = size(Y,2);
        N{i} = size(Y,1);
        indTr = globalOpt.indTr;
        if indTr == -1
            indTr = 1:N{i};
        end
        if ~exist('Yts')
            indTs = setdiff(1:size(Y,1), indTr);
            Yts{i} = Y(indTs,:);
        end
        Ytr{i} = Y(indTr,:);

        t{i} = linspace(0, 2*pi, size(Y, 1)+1)'; t{i} = t{i}(1:end-1, 1);
        timeStampsTraining{i} = t{i}(indTr,1); %timeStampsTest = t(indTs,1);
    end

    for i=2:numberOfDatasets
        if N{i} ~= N{i-1}
            error('The number of observations in each dataset must be the same!');
        end
    end

    %%
   % if ~onlyTest
        % Free up some memory
        clear('Y')

        options = svargplvmOptions(Ytr, globalOpt, labelsTrain);



        if ~isempty(globalOpt.dynamicsConstrainType)
            %for i=1:numberOfDatasets
                % Set up dynamcis (i.e. back-constraints) model
                optionsDyn.type = 'vargpTime';
                optionsDyn.inverseWidth=30;
                %   optionsDyn.vardistCovars = vardistCovarsMult;
                optionsDyn.initX = globalOpt.initX;
                optionsDyn.constrainType = globalOpt.dynamicsConstrainType;

                if exist('timeStampsTraining')
                    optionsDyn.t = timeStampsTraining;
                end
                if exist('labelsTrain') && ~isempty(labelsTrain)
                    optionsDyn.labels = labelsTrain;
                end
            %end
        else
            optionsDyn= [];
        end




        model = svargplvmModelCreate(Ytr, globalOpt, options, optionsDyn);
        if exist('diaryFile')
            model.diaryFile = diaryFile;
        end

        model.globalOpt = globalOpt;
        model.options = options;

        %%%% TEMP
        if exist('whiteVar')
            model.dynamics.kern.comp{2}.variance = whiteVar;
        end
        %%%%

        % Force kernel computations
        params = svargplvmExtractParam(model);
        model = svargplvmExpandParam(model, params);

        %%

        if displayIters
            model = svargplvmOptimiseModel(model);
        else
            model = svargplvmOptimiseModelNoDisplay(model);
        end
        
    %end

    %%

    capName = model.globalOpt.dataSetName;
    capName(1) = upper(capName(1));
    errors = fgplvmNearestNeighbour(model.comp{1},lblsTrain);
    model2 = vargplvmReduceModel(model.comp{1}, 2);
    errors2 = fgplvmNearestNeighbour(model2, lblsTrain);
    fprintf('# Visualisation errors in all dims/in 2D:  %d / %d\n', errors, errors2)
    if printPlots
        vargplvmPrintPlot(model2, lblsTrain, [capName 'Vargplvm'], model.globalOpt.experimentNo);

        %%
        % For 3D plots:
        % labels = model.options{1}.labels; dims = [1 3 5];
        % plot3k({model.X(:,dims(1)) model.X(:,dims(2)) model.X(:,dims(3))}, 'ColorData', labels, 'Marker', {'x',6});
    end

    %%
    obsMod = 1; % one of the involved sub-models (the one for which we have the data)
    infMod = setdiff(1:2, obsMod);

    % sharedDims = svargplvmFindSharedDims(model);

    %---------------------------- PREDICTIONS ---------------
    if ~doPredictions
        return
    end


    %%
    fprintf('\n# PREDICTIONS: \n\n');
    if ~exist('testOnTraining')
        testOnTraining=0;
    end


    [sharedDims, privateDims] = svargplvmFindSharedDims(model,[],[],{obsMod infMod});

    
    % svargplvmPredictions
    [x_star_all, varx_star_all, mini] = vargplvmPredictLatent(model.comp{obsMod}, Yts{obsMod}, [], false, model.globalOpt.reconstrIters,0,[],[],1);
    infMethod = 1; % Experiment with different inference methods as well (check svargplvmPredictionsFunc)
    [Zpred, testInd, XpredAll, varXpredAll, indNN] = ...
        svargplvmPredictionsFunc(model, 0, x_star_all, varx_star_all, obsMod, infMod, [], 1, infMethod);
    % Store into a matrix
    ZpredMuAll = zeros(length(Zpred), size(Zpred{1},2));
    for i = 1:length(Zpred)
        ZpredMuAll(i,:) = Zpred{i};
    end

    %--- For classification we are interested in labels
    ZpredMuAll(ZpredMuAll >= 0) = 1;
    ZpredMuAll(ZpredMuAll < 0) = -1;
    % This dataset should have exactly one label present, so if all zeros
    % are returned take the maximum ind. as one (this in practise never
    % happens)
    tmpList = find(sum(ZpredMuAll,2) == -3);
    for tt = 1:length(tmpList)
        [~,tmpInd] = max(Zpred{tmpList(tt)});
        ZpredMuAll(tmpList(tt),tmpInd) = 1;
    end
    
    NNpred = Ytr{infMod}(mini,:);
    realLabels = lblsTest(testInd,:);

    % Find the error in the labels ------------
    [labelErrors, labelSharingErrors, wrongIndices] = findLabelErrors(realLabels, ZpredMuAll);
    [labelErrorsNN, void, wrongIndicesNN] = findLabelErrors(realLabels,NNpred);

    fprintf('# Gplvm label errors: %d\n', labelErrors);
    fprintf('# Gplvm label sharing errors: %d\n', labelSharingErrors);
    fprintf('# NN label errors:    %d\n', labelErrorsNN);

    try
        % Confusion matrix
        cmat_a = confMatrix( XpredAll, model.X, labelsTrain, labelsTest, 1);
        cmat_rel = cmat_a ./ repmat(sum(cmat_a, 2), 1, length(unique(labelsTrain)));
       % cmat_rel
        cmat_relAll = cmat_relAll + cmat_a;
        labelErrorsAll = labelErrorsAll + labelErrors;
        labelErrorsNNAll = labelErrorsNNAll + labelErrorsNN;
    catch e
        getReport(e)
    end

    errors.gplvm = labelErrorsAll;
    errors.NN = labelErrorsNNAll;
    errors.wrongIndicesGplvm = wrongIndices;
    errors.wrongIndicesNN = wrongIndicesNN;
    errors.confusionMatrix = cmat_relAll;
    predictions.Z = ZpredMuAll;
    predictions.X = XpredAll;
    predictions.varX = varXpredAll;
    scales = svargplvmShowScales(model,0);
    save(['demOilSvargplvmPred' num2str(experimentNo) '.mat'], 'errors', 'predictions','scales');
    
%     
%     
%     
%     scales1All{trialNo} = model.comp{1}.kern.inputScales;
%     scales2All{trialNo} = model.comp{2}.kern.inputScales;
%     
%     
%     indsAll{trialNo} = indTrial;
%     labelErrorsAll{trialNo} = labelErrors;
%     labelSharingErrorsAll{trialNo} = labelSharingErrors;
%     labelErrorsNNAll{trialNo} = labelErrorsNN;
%     wrongIndicesAll{trialNo} = wrongIndices;
%     wrongIndicesNNAll{trialNo} = wrongIndicesNN;
%     SNR1{trialNo} = (1/model.comp{1}.beta)/var(model.comp{1}.m(:));
%     SNR2{trialNo} = (1/model.comp{2}.beta)/var(model.comp{2}.m(:));
%     cmat_relAll{trialNo} = cmat_rel;
    
    
   % % open the file with write permission
   % fid = fopen(['demOilSvargplvmMultiple' num2str(experimentNo) '.txt'], 'w+');
   % fprintf(fid, 'scales1: );
   % fclose(fid);
    
    %------------------------------------------------------------------
    
    % Temporary ugly solution... when the experiments are over just forget
    % about the "allDataPerClass'.
    if  exist('allDataToKeep')
        flag = 1;
        allDataPerClass = allDataToKeep;
    else
        flag = 0;
    end
    keep('flag','experimentNo', 'jj','allDataPerClass',...
        'curIndAll','globalOpt','doPredictions','printPlots','labelsTest', ...
        'YtrFull','lblsTrainFull','labelsTrainFull','Yts','lblsTest','trialNo','totalTrials', ...
        'labelErrorsAll', 'labelErrorsNNAll','cmat_relAll');
    
    if flag
        allDataToKeep = allDataPerClass;
    end
end










