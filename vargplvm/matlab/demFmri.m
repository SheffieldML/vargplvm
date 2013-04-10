enableParallelism = 0;
delete LOG_fmri400Processed3_1.txt
diary LOG_fmri400Processed3_1.txt
try
    load fmri400Processed3
    experimentNo = 1;%globalExperimentNumber('vargplvm','fmri','get');
    indPoints = 120;
    latentDim = 60;
    dataSetSplit = 'custom';
    indTr = 1:size(Y,1);
    indTs = size(Y,1);
    dynUsed = 0;
    initVardistIters = 250;
    dataSetName = 'fmri4003_';
    itNo = [1000 1000 4000 500 500 500 500];
    tic
    demHighDimVargplvm3
    a=toc
    diary off
catch
    diary off
    toc
end


% Note: By only doing the initialisation, I get 3 scales.