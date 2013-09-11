% This demo compares dynamic Bayesian GPLVM with and without image pyramids
% to see if we get a better reconstruction of the dog tail motion using
% pyramids.
%
% Created by Andreas and Teo on 31 July 2013, based on demosDynamics.m

%------------ GENERATION -------
% Train

%% Setting paths
addpath(genpath('../../../GPmat/matlab/')); % GPmat
% addpath(genpath('../../../vargplvm/vargplvm/matlab/')); % Bayesian GPLVM
addpath('../../../GPmatlab/netlab3_3/'); % NetLab, which is not included in gitHub

%% Preparing pyramid dataset:
% This only needs to be done once.
load dog
% To view it:
% for f=1:600, image(reshape(Y(mod(f-1,61)+1,:),[360 640])); colormap gray(256); pause(0.033); end
levels = 5;
Y = im2pyramid(Y, height, width, levels);
save('dogPyramid', 'Y', 'height', 'width', 'levels');
clear 

%% Selecting dataset -- To run experiments, change the dataset below and run the test again:
% function doAll(dataSetName)
dataSetName = 'dog';
% Running the rest for each dataset. 
fprintf(1,'\n\n#-----  DOG DEMO: Generation ----#\n');
experimentNo=61;
indPoints=-1; 
latentDim=35;
fixedBetaIters= 0;
initVardistIters = 800;
mappingKern = 'rbfardjit';
reconstrIters = 1; % no reconstruction needed here
itNo=repmat(1000, 1, 16); %16000
periodicPeriod = 4.3983; % Recalculated for dataToKeep=60
dynamicKern={'rbfperiodic','whitefixed','bias','rbf'};
vardistCovarsMult=0.8;
whiteVar = 1e-6;
dataToKeep = 60; 
dataSetSplit = 'custom';
indTr = [1:60];
indTs = 60; % We don't really reconstruct in this experiment
learnSecondVariance = 0;
DgtN = 1; % Uncomment for the faster version
predWithMs = 0;
%%
demHighDimVargplvm3

%%
%-- Then generate:
% experimentNo=61; 
% dataToKeep = 60; 
% dataSetSplit = 'custom';
% indTr = [1:60]; 
% indTs = 60;
futurePred = 40; 
doSampling = 0; 
demHighDimVargplvmTrained
%clear Yts; clear YtsOriginal; clear Testmeans2; clear Testcovars2;
%playMov(height, width, [], [Ytr(end-5:end,:); Varmu2]);

%%

% Produce plots
bar(prunedModelInit.kern.comp{1}.inputScales)
print -depsc ../diagrams/dog_scalesInit.eps; system('epstopdf ../diagrams/dog_scalesInit.eps');
bar(model.kern.comp{1}.inputScales)
print -depsc ../diagrams/dog_scalesOpt.eps; system('epstopdf ../diagrams/dog_scalesOpt.eps');

nPixels = height*width; % This will be used in case Y was built using pyramids.
fr=reshape(Ytr(end,1:nPixels),height,width); imagesc(fr); colormap('gray'); % Last training image
print -depsc ../diagrams/dogGeneration_lastOfTraining.eps; system('epstopdf ../diagrams/dogGeneration_lastOfTraining.eps');
fr=reshape(Varmu2(1,1:nPixels),height,width); imagesc(fr); colormap('gray');  % First predicted
print -depsc ../diagrams/dogGeneration_firstOfTest.eps; system('epstopdf ../diagrams/dogGeneration_firstOfTest.eps');
fr=reshape(Varmu2(13,1:nPixels),height,width); imagesc(fr); colormap('gray'); % A subsequent frame
print -depsc ../diagrams/dogGeneration_frame14.eps; system('epstopdf ../diagrams/dogGeneration_frame14.eps');


for i=1:size(Varmu2,1)
    img = reshape(Varmu2(i,1:nPixels), height, width);
    % imagesc(img); colormap('gray'); 
    imwrite(uint8(img), sprintf('../diagrams/%s_%02d.png',dataSetName, i))
    % pause(0.1);
end

%% Compare the dogs
diffs = zeros(40,2);
for f=1:40,
    imgFlat = imread(sprintf('../diagrams/dog_%02d.png',f));
    imgPyramid = imread(sprintf('../diagrams/dogPyramid_%02d.png',f));
    diffImg = abs(imgFlat - imgPyramid);
    imagesc(diffImg); 
    colorbar
    diffs(f,1) = mean(diffImg(:));
    diffs(f,2) = max(diffImg(:));
    fprintf('f=%d, mean_difference=%f, max_difference=%f\n', f, diffs(f,1), diffs(f,2));
    pause(0.1);
end
