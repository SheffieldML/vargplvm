
%###------------   Predictions script ---------------------------------####
allExperiments = [201, 5554]; % lovelace: 5-7, compbio4: 9,11
warning off
for experimentNo = allExperiments
	try
		%diaryFile=['LOGS/LOG_demMultiLabelEmotions' num2str(experimentNo) '.txt']; diary(diaryFile);
		dataSetName = 'mulan_emotions';
		onlyTest = 1;
        doPredictions = 1;
		demClassificationGeneral
       % save(['demMulan_emotionsSvargplvmPred' num2str(globalOpt.experimentNo)],'ZpredMuAll')
	catch e
	    warning(['Problem with experiment' num2str(experimentNo)])
		getReport(e)
	end
	%diary off
	keep('experimentNo','allExperiments')
end
warning on


%###------------   Eval. script ---------------------------------####
allExperiments = [18 19 20]; % lovelace: 5-7, compbio4: 9,11
for experimentNo = allExperiments
	try
		dataSetName = 'mulan_emotions';
		onlyTest = 1;
		demClassificationGeneral
		doPredictions = false;
		svargplvmShowScales(model,false)
		svargplvmShowScales(model);
		title(num2str(experimentNo))
		sharedDims = svargplvmFindSharedDims(model)
		for j=1:length(model.comp)
            fprintf('# SNR%d = %.6f\n',j,(1/model.comp{j}.beta)/var(model.comp{j}.m(:)));
        end
		pause
	catch e
	    warning(['Problem with experiment' num2str(experimentNo)])
		getReport(e)
	end
	keep('experimentNo','allExperiments')
end



%###------------   Reoptimise script ---------------------------------####
allExperiments = [18]; % lovelace: 5-7, compbio4: 9,11
for experimentNo = allExperiments
	try
		diaryFile=['LOGS/LOG_demMultiLabelEmotions' num2str(experimentNo) '.txt']; diary(diaryFile);
		dataSetName = 'mulan_emotions';
		onlyTest = 1;
		demClassificationGeneral
		model = svargplvmReOptimiseModel(model, 0, 500, experimentNo*100);
		svargplvmShowScales(model,false)
		svargplvmShowScales(model);
		sharedDims = svargplvmFindSharedDims(model)
	catch e
	    warning(['Problem with experiment' num2str(experimentNo)])
		getReport(e)
	end
	diary off
	keep('experimentNo','allExperiments')
end