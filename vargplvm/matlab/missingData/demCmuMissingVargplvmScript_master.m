clear;

ss = 2; % seed
missingPerc = 0.1:0.05:0.5;
missingPerc = [missingPerc 0.6:0.1:0.9];
experimentNo = 1;
fprintf('## EXPERIMENT NO = %d\n\n\n', experimentNo)

errorMeanAll = nan(1, length(missingPerc));
errorBGPLVM_missingAll = nan(1, length(missingPerc));
           
model = [];

for mm = 1:length(missingPerc)
    keep('model', 'experimentNo','ss','mm', 'missingPerc','errorMeanAll', 'errorBGPLVM_missingAll');
    rseed = ss;
    QmissingPerc = missingPerc(mm);
    
    Ntr = 50;
    Ninterp = 200;
    Nmissing = 80;
    dataSetName = 'cmu35gplvm';
    numActive = 35;
    iters = [1500 1000];
    initIters = 3000;
    initSNR = 80;
    init_missing = 'posterior';
    linearMapping = false;
    
    displayIters = false;
    runBGPLVMmisFixed = false;
    runGP = false;
    if mm == 1
        runBGPLVM = true;
    else
        runBGPLVM = false;
    end
    demRegressionVargplvm2_master
    
    errorMean = mean(abs(Yinterp(:) - meanPred(:)));
    errorBGPLVM_missing = mean(abs(Yinterp(:) - YpredBGPLVM_missing(:)));
    
    
    errorMeanAll(mm) = errorMean;
    errorBGPLVM_missingAll(mm) = errorBGPLVM_missing;
    
end


return



%% PLOT

h=figure;
lw = 1.5;
ms = 8;
fs = 13;

errorMeanAll = errorMeanAll(ind,xind);
errorBGPLVM_missingAll = errorBGPLVM_missingAll;
plot(missingPerc(xind), mean(errorBGPLVM_missingAll,1), ['kd' '-'], 'LineWidth', lw, 'MarkerSize', ms);

