%cd ../../vargplvmMex/
%load tmpData
%delete DEBUG.txt

writeMatrixFile('means.txt', vardist.means);
writeMatrixFile('covars.txt',vardist.covars);
writeMatrixFile('asPlus1.txt',asPlus1);
writeMatrixFile('aDasPlus1.txt', aDasPlus1);
writeMatrixFile('ZmZm.txt', ZmZm);
writeMatrixFile('covGrad.txt', covGrad);

ZmZm_arr = size(ZmZm,2);
covGrad_arr = size(covGrad,2);
% mex vargplvm.cpp
vargplvm([M,N,Q],A,ZmZm_arr, covGrad_arr)


partInd2= dlmread('partInd2.txt',  '\t');
partInd2 = partInd2';
partA2= dlmread('partA2.txt',  '\t');
partA2 = partA2';
gVarmeans= dlmread('gVarmeans.txt',  '\t');
gVarmeans = gVarmeans';
gVarcovars= dlmread('gVarcovars.txt', '\t');
gVarcovars = gVarcovars';
        