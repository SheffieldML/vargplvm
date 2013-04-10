% DEMCMU35VARGPLVMRECONSTRUCTTAYLOR Reconstruct right leg and body of CMU 35.
%
% DESC
%  Load a model learned with vargplvm and reconstruct the missing indices
%  using the Taylor framework developed for fgplvm (angle errors)
% Based on demCmu35fgplvmReconstructTaylor.m

% See also: vargplvmTaylorAngleErrors
% VARGPLVM

randn('seed', 1e5);
rand('seed', 1e5);


% dynUsed=1;
%
% % Get the sequence numbers.
 [Y, lbls] = lvmLoadData('cmu35WalkJog'); %%%
 seq = cumsum(sum(lbls)) - [1:31]; %%%

dataSetName = 'cmu35gplvm';

if ~exist('experimentNo') experimentNo = 1; end

% load data
[Y, lbls, Ytest, lblstest] = lvmLoadData(dataSetName);



dataSetName = 'cmu35gplvm';

capName = dataSetName;
capName(1) = upper(capName(1));

fileName=['dem' capName 'Vargplvm' num2str(experimentNo)];

% load data
[Y, lbls, Ytest, lblstest] = lvmLoadData(dataSetName);

%---- Remove some sequences %%%%%%
seqFrom=7;
seqEnd=7;


if seqFrom ~= 1
    Yfrom = seq(seqFrom-1)+1;
else
    Yfrom = 1;
end
Yend=seq(seqEnd);
Y=Y(Yfrom:Yend,:);
seq=seq(seqFrom:seqEnd);
seq=seq - ones(1,length(seq)).*(Yfrom-1);
%-----

origBias = mean(Y);
origScale = 1./sqrt(var(Y));
% %scale = ones(size(scale));
% Y = Y - repmat(origBias, size(Y, 1), 1);
% Ytest = Ytest - repmat(origBias, size(Ytest, 1), 1);
% Y = Y.*repmat(origScale, size(Y, 1), 1);
% Ytest = Ytest.*repmat(origScale, size(Ytest, 1), 1);

% REMOVE LEG TEST
% Indices associated with right leg.
legInd = [8:14];
% Indices associated with upper body.
bodyInd = [21:50];


dt=0.05;
timeStampsTest = ([0:size(Ytest,1)-1].*dt)';


%for experimentNo = 1:3;

% Load saved model.
%load(['dem' capName num2str(experimentNo) '.mat']);
load(fileName); % From demCmu35gplvmVargplvm1.m
fprintf(1,'# Loaded %s\n',fileName);
if isfield(model, 'dynamics') & ~isempty(model.dynamics)
    fprintf(1,'# Dynamics kernel:');
    kernDisplay(model.dynamics.kern);
    model.dynamics.t_star = timeStampsTest;
    kernName = model.dynamics.kern.comp{1}.type;
else
    kernName = 'static';
end
fprintf(1,'# experimentNo=%d\n',experimentNo);
fprintf(1,'# Reconstruction Iters=%d\n',model.reconstrIters);

 model.reconstrIters = 2000; %%%TEMP



startInd = 63;

legErrs = vargplvmTaylorAngleErrors(model, Y, Ytest, startInd, origBias, origScale, legInd, [dataSetName ...
    'Leg'], experimentNo);

save(['dem' capName 'ReconstructLegs' num2str(experimentNo) '.mat'], 'legErrs');

%----
figPath = ['Results/CMU/' num2str(experimentNo) kernName '/Leg/'];
%try
%    saveAllOpenFigures(figPath, [],1)
%catch 
%    %
%end
%----

legErrs



bodyErrs = vargplvmTaylorAngleErrors(model, Y, Ytest, startInd, origBias, origScale, bodyInd, [dataSetName ...
    'Body'], experimentNo);

save(['dem' capName 'ReconstructBody' num2str(experimentNo) '.mat'], 'bodyErrs');


figPath = ['Results/CMU/' num2str(experimentNo) kernName '/Body/'];
%try
%    saveAllOpenFigures(figPath, [],1)
%catch 
%    %
%end

%----
bodyErrs

% bar(model.kern.comp{1}.inputScales)


%end


