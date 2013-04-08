distcomp.feature( 'LocalUseMpiexec', false )
matlabpool
try
load fmriProcessed900
experimentNo = 2;%globalExperimentNumber('vargplvm','fmri','get');
indPoints = 115;
dataSetSplit = 'custom';
indTr = 1:336;
indTs = 336;
dynUsed = 0;
initVardistIters = 300;
dataSetName = 'fmri900';
itNo = [1000 500];
demHighDimVargplvm3
diary off
exit
catch e
    diary off
    exit
end