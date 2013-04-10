function printGlobalOptions(globalOpt)

% PRINTGLOBALOPTIONS
% VARGPLVM

fprintf(1,'\n#----------------------------------------------------\n');
fprintf(1,'# Dataset: %s\n',globalOpt.dataSetName);
fprintf(1,'# ExperimentNo: %d\n', globalOpt.experimentNo);
fprintf(1,'# Latent dimensions: %d\n',globalOpt.latentDim);
fprintf(1,'# Iterations (with/without fixed Beta): %d / %s\n',globalOpt.fixedBetaIters, num2str(globalOpt.itNo));
fprintf(1,'# Reconstruction iterations: %d\n',globalOpt.reconstrIters);
fprintf(1,'# Tie Inducing points: %d\n',globalOpt.fixInd);
fprintf(1,'# InitX: %s\n',globalOpt.initX);
fprintf(1,'# Reoptimise inducing points (for reconstr): %d \n',globalOpt.testReoptimise);
fprintf(1,'# Dynamics used: %d\n', globalOpt.dynUsed);
if dynUsed
    fprintf(1,'# Dynamics kern: ');
    disp(globalOpt.dynamicKern);
end

fprintf(1,'# VardistCovarsMult: %d \n', globalOpt.vardistCovarsMult);
fprintf(1,'# InvWidthMultDyn: %d \n', globalOpt.invWidthMultDyn);
fprintf(1,'# InvWidthMult: %d \n', globalOpt.invWidthMult);
printf(1,'#----------------------------------------------------\n');
