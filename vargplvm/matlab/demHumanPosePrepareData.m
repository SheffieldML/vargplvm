% DEMHUMANPOSEPREPAREDATA Load the human pose dataset and prepare the final dataset
% by subsampling or selecting a subset of the sequences.
% DESC Load the human pose dataset and prepare the final dataset
% by subsampling or selecting a subset of the sequences.
%
% COPYRIGHT: Andreas C. Damianou, Carl Henrik Ek,  2011
% SEEALSO : demHumanPoseSvargplvm1
%
% VARGPLVM



% Sequences are:
% 1) Walk   -->
% 2) Drunk
% 3) Walk -->
% 4) Walk <--
% 5) Walk in circle
% 6) Walk in circle
% 7) Something weird: going a bit backwards and then standing
% 8) Walk in an ellipsis

%{
load sil_images
Yim = Y;
Yim_test = Y_test;

load humanPose;

%ind2 = floor(1:2:size(Y,1));
%%%__silhouettes are whole images. If imageSilhouette is true, then all
%%%pixels will be used (instead of just features extracted).
if exist('imageSilhouette') && imageSilhouette
    fprintf('# Using the actual silhouette''s pixels!\n');
    Y = Yim;
    Y_test = Yim_test;
    %clear Yim
    %clear Yim_test
end
%{
% Remove the 'drunk' sequence (2nd seq).
Y1 = Y(1:100,:);
Y2 = Y(337:end,:);
Y = [Y1; Y2];
clear('Y1','Y2');
Z1 = Z(1:100,:);
Z2 = Z(337:end,:);
Z = [Z1; Z2];
clear('Z1','Z2');
%}

% Subsample
Yall{1} = Y;
Yall{2} = Z;


%--- Uncomment this code to get a subset of the sequences
%  seqFrom=3; %3,4
%  seqEnd=4;
% 
%  if seqFrom ~= 1
%      Yfrom = seq(seqFrom-1)+1;
%  else
%      Yfrom = 1;
%  end
%  Yend=seq(seqEnd);
%  Yall{1}=Yall{1}(Yfrom:Yend,:);
%  Yall{2}=Yall{2}(Yfrom:Yend,:);
%  seq=seq(seqFrom:seqEnd);
%  seq=seq - ones(1,length(seq)).*(Yfrom-1);
%---
subSamp = 1;
if exist('inds')
    seqOrig = seq;
    %inds = [{1:100} {337:405} {406:471} {955:1015}];
    seq = []; indsAll = []; prev = 0;
    for i=1:length(inds)
        inds{i} = inds{i}(1:subSamp:end);
        indsAll = [indsAll inds{i}];
        prev = prev + length(inds{i});
        seq = [seq prev];
    end
    Yall{1} = Yall{1}(indsAll,:);
    Yall{2} = Yall{2}(indsAll,:);
    Yim = Yim(indsAll,:);
end
%}

%--------------------------------------------------------

[Yall, lbls, YtestAll] = svargplvmLoadData('humanPoseAll');
Y = Yall{1};
Yim = Yall{2};
Z = Yall{3};

if ~exist('testSeq')
    Y_test = YtestAll{1};
    Yim_test = YtestAll{2};
    Z_test = YtestAll{3};
else
    Y_test = Yall{1}{testSeq};
    Yim_test = Yall{2}{testSeq};
    Z_test = Yall{3}{testSeq};
end




height = lbls{1}; width = lbls{2}; seqOrig = lbls{3};


% Subsample sequences
seq = []; prev = 0;

if ~exist('subsample')
    subsample = 3;
end
for i=1:length(Y)
    Y{i} = Y{i}(1:subsample:end,:);
    Yim{i} = Yim{i}(1:subsample:end,:);
    Z{i} = Z{i}(1:subsample:end,:);
    tempN = size(Y{i},1);
    seq = [seq tempN+prev];
    prev = prev+tempN;
end

% Keep only a few sequences
Yall2 = [];
YimAll = [];
Zall = [];

if ~exist('seqToKeep')
    seqToKeep = [1 3 4 5 6 7 8];
end
seq = []; prev = 0;
for i=1:length(seqToKeep)
    Yall2 = [Yall2;Y{seqToKeep(i)}];
    YimAll = [YimAll;Yim{seqToKeep(i)}];
    Zall = [Zall;Z{seqToKeep(i)}];
    tempN = size(Y{seqToKeep(i)},1);
    seq = [seq tempN+prev];
    prev = prev+tempN;
end
Yim = YimAll;
clear Yall


if exist('usePixels') & usePixels
    Yall{1} = Yim;
    Yall_tr1 = Yall2;
else
    Yall{1} = Yall2;
end
Yall{2} = Zall;

clear('YimAll','Yall2','Zall','Y','Z','YtestAll','prev');


% Make sure that sequences are correct
%{
for j=1:length(seq)
    k=1;
    for i=seq(j)-3:seq(j)+3
        subplot(2,4,k)
        imagesc(reshape(Yim(i,:), height, width))
        k = k+1;
    end
    pause
end
%}


% Print all sequences
%{
h1=subplot(1,2,1);
handle=xyzankurVisualise(Z_test(1,:), 1);
h2=subplot(1,2,2);
imagesc(reshape(Yim_test(1,:),128,128))
colormap('gray')
for i=2:size(Yim_test,1)
    h1=subplot(1,2,1);
    %xyzankurModify(handle, Z_test(i,:));
    xyzankurVisualise(Z_test(i,:), 1);

    h2=subplot(1,2,2);
    imagesc(reshape(Yim_test(i,:),128,128))
    pause(0.1)
    clf
end
%}