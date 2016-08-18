seeds = 2:5;
missingPerc = 0.05:0.05:0.5;
missingPerc = [missingPerc 0.6:0.1:0.9];
if ~exist('experimentNo','var'), experimentNo = 106; end


errorMeanAll = nan(length(seeds), length(missingPerc));
errorBGPLVM_missingAll = nan(length(seeds), length(missingPerc));

            
model = [];
for ss = 1:length(seeds)
    for mm = 1:length(missingPerc)
        keep('model', 'experimentNo','ss','seeds','mm', 'missingPerc','errorLinAll','errorNNAll', 'errorNNxAll','errorMeanAll', 'errorBGPLVMAll', 'errorBGPLVM_missingAll', 'errorBGPLVM_fixed_missingALL', 'errorGPAll', 'errorGPLVMAll','errorBGPLVMInitAll');   
        rseed = seeds(ss);
        QmissingPerc = missingPerc(mm);  
               
        Ntr = 40;
        D = 5;
        Q = 15;
        Ninterp = 100;
        Nmissing = 60;
        outputNoise = 0.001;
        numActive = 30;
        iters = 2000;
        initIters = 2000;
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
        runGPLVM = false;
        
        demRegressionVargplvm2_master
        
        errorMean = mean(abs(Yinterp(:) - meanPred(:)));
        errorBGPLVM_missing = mean(abs(Yinterp(:) - YpredBGPLVM_missing(:)));        
        errorMeanAll(ss,mm) = errorMean;
        errorBGPLVM_missingAll(ss, mm) = errorBGPLVM_missing;
    end
end


