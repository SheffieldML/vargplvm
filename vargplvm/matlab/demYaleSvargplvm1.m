% DEMYALESVARGPLVM1 Run the Shared Var. GP-LVM on a subset of the Yale
% faces.
% DESC Run the Shared Var. GP-LVM on a subset of the Yale faces. The code
% for creating this subset out of raw images exists in comments. This demo
% us a wrapper for the generic shared var. GP-LVM demo
% (demSharedVargplvm1). 

% SVARGPLVM

% Create the dataset out of the images.
baseDir=[localDatasetsDirectoryLarge 'CroppedYale' filesep 'CroppedYale'];
selDirs = {'04','07','26','31', '19','30'};

for d=1:length(selDirs)
    dirFrom=[baseDir filesep 'yaleB' selDirs{d}];
    a=dir(dirFrom);
    counter = 0;
    for i=1:length(a)
        if length(a(i).name)>4 & strcmp(a(i).name(end-2:end),'pgm') ...
                & ~strcmp(a(i).name(end-10:end-4),'Ambient')
            im = imread([dirFrom filesep a(i).name]);
            %imagesc(im), colormap('gray'); title(a(i).name), pause
            counter = counter+1;
            Yall{d}(counter,:)=im(:)';
        end
    end
    Yall{d} = double(Yall{d});
end
height = size(im,1);
width = size(im,2);
dataSetsUsed = selDirs;
numberOfDatasets = length(Yall);
Y = Yall;

for i=1:size(Y{1},1)
    for d=1:numberOfDatasets
        subplot(1,numberOfDatasets,d)
        imagesc(reshape(Y{d}(i,:),height, width)), colormap('gray');
    end
    pause
end   
keep('Y','dataSetsUsed','height','width');

return
%%


experimentNo = 1;
dataType = 'Yale4Sets';
dataSetNames = 'YaleSubset4_2';
enableParallelism = 0;
[Y,lbls]=svargplvmLoadData(dataSetNames);
Yall{1} = [Y{1};Y{2}];
Yall{2} = [Y{3};Y{4}];
clear Y;

numberOfDatasets = length(Yall);
height = lbls(1); width = lbls(2);

%{
for i=1:size(Yall{1},1)
    for d=1:numberOfDatasets
        subplot(1,numberOfDatasets,d)
        imagesc(reshape(Yall{d}(i,:),height, width)), colormap('gray');
    end
    pause
end   
%}


% TODO
initLatent = 'ppca';
indPoints = 100;
latentDimPerModel = 10;
%mappingKern = 'rbfardjit';
initVardistIters = 100;
itNo = [1000 1000 1000 500 500 500 500] ;

demSharedVargplvm1
