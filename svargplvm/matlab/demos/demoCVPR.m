% DEMOCVPR Commands to recreate the experiments for the CVPR paper.
% DESC This is a 'hyper-demo', that calls other demos.
% The following commands set some hyper-parameters of the algorithm and
% call the corresponding demo to run the experiment.
% In particular, 'itNo' corresponds to the number of iterations for the
% optimiser. We have observed that even a much smaller number than the default
% does not change the results noticeably.
%
% COPYRIGHT: Andreas C. Damianou, Carl Henrik Ek,  2011
%
% SHEFFIELDML


%--

%----- Yale faces:

% a) Dataset
% The called demo also creates the dataset out of the image files. You only need to specify
% the path to the folder containing the images.
%  We use ONLY the images of the 
% basic pose i.e. all images should start with: yaleBXX_P00A-, where XX
% is the number of the subject (see demo for which subjects are used).

% b) Demo
latentDimPerModel = 7;
itNo = [1000 1000 1000 1000 1000]; 
mappingKern = 'rbfardjit';
DgtN = 1;
initial_X = 'together';
indTr = [1:192]
demYaleSvargplvm2
% The pre-trained model is experimentNo = 25

%---- Human Pose:
% a) Dataset
% The dataset used is the one mentioned in the paper but here we only take
% a subset. The original dataset that needs to be provided contains 8
% sequences, having 100, 236, 69, 66, 483, 577, 238 and 158 frames
% respectively. The demo preprocesses the data (i.e. subsamples every 2
% frames since the data are very dense and only keeps a subset of the
% sequences) assuming that they come in the format described in the
% beginning of 'demHumanPosePrepareData.m'. In short, a .mat file needs to
% be placed in the path, and this .mat file should contain 3 cells, each
% corresponding to the HoG, pixel and pose features respectively. Each of
% these representations is itself a 1x8 cell, i.e. one cell with the
% corresponding features for each sequence.
% The demo discards sequences 2 and 7, because they correspond to motions
% very dissimilar to the rest and somehow unnatural. The 8th sequence is
% used as a test set. We don't use the provided by the authors test set, as
% it is very very similar to one of the training sequences.


% b) Demos



% Human pose, dynamical model:
experimentNo = 120; % identical: 10000, but 120 is the original run
subsample = 2; % Experiment id: #120 (experiment id #10000 is identical as well)
seqToKeep = [1 3 4 5 6];
testSeq = 8;
initial_X = 'separately';
latentDimPerModel = [9 6];
dynUsed=1;
initVardistIters = 240;
itNo = [800 700 800 800];
vardistCovarsMult = 0.09;
indPoints = 100;
doPredictions = 0; 
demHumanPoseSvargplvm1
clear
experimentNo = 120; dynUsed=1; demHumanPosePredictions 
%(Something seems to be wrong with the current implementation,it finds 1 shared dim, not 2 as in 120)


% Human pose, non-dynamical model:

subsample = 2; % Experiment id: #311/217
seqToKeep = [1 3 4 5 6];
testSeq = 8;
initial_X = 'separately';
initVardistIters = 240;
mappingKern = 'rbfardjit';
latentDimPerModel = [5 3];
dynUsed=0;
itNo = [400 400 400 500 500];
indPoints = 110;
demHumanPoseSvargplvm1
