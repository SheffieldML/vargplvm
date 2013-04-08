% Before using that set vargplvmOptimise to do gradChek and also either
% comment the part where tha variance for theta_t is held fixed or be
% careful to see that only in that element the gradient is incorrect.

clear; itNo=1; indPoints=15; latentDim=5; fixInd=1; fixedBetaIters=0; dynUsed=1; demHighDimVargplvm1
load delta; load gradient; load deltaf;
d=norm(deltaf - gradient)/norm(gradient + deltaf); %%
d1=norm(deltaf - gradient,1)/norm(gradient + deltaf,1); %%
fprintf(1,' Norm1 difference: %d\n Norm2 difference: %d\n',d1,d);

