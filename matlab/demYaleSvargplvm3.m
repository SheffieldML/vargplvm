
% TODO!!!! 
% This demo is like demYaleSvargplvm2, but it adds a few pairs in the end,
% where the first part of the pair is a NEW face, and the second part is
% one of the training ones but with the same angle. The NEW face does NOT
% come with every possible angle, so that we have to fill it in after
% training via generation (latent space sampling)


latentDimPerModel = 7;
itNo = [1000 1000 1000 1000 1000]; 
mappingKern = 'rbfardjit';
DgtN = 1;
initial_X = 'together';
indTr = [1:192];


%% < LOAD DATA FROM demYaleSvargplvm2 >


% Create the dataset out of the images.
baseDir=[localDatasetsDirectoryLarge 'CroppedYale' filesep 'CroppedYale'];
selDirs = {'37'};

    dirFrom=[baseDir filesep 'yaleB' selDirs{1}];
    a=dir(dirFrom);
    counter = 0;
    for i=1:length(a)
        if length(a(i).name)>4 & strcmp(a(i).name(end-2:end),'pgm') ...
                & ~strcmp(a(i).name(end-10:end-4),'Ambient')
            im = imread([dirFrom filesep a(i).name]);
            %imagesc(im), colormap('gray'); title(a(i).name), pause
            counter = counter+1;
            Ynew(counter,:)=im(:)';
        end
    end
    Ynew = double(Ynew);

height = size(im,1);
width = size(im,2);

indsToSelect = [1:10];
Yextra{1} = [];
Yextra{2} = [];
for i=1:length(indsToSelect)
    % The seecond dataset has as extra face the new face
    Yextra{2} = [Yextra{2} ; Ynew(i,:)];
   % Yextra{2} = [Yextra{2} ; Ytr{2}(i,:)];
    
    % For the first dataset, we repeat any consistent face from the rest of
    % the Ytr, any modality.
    % Choose modality
    m = randperm(numberOfDatasets);
    m = m(1);
    % Choose person within modality
    person = randperm(length(unique(identities{m})));
    person = person(1);
    allAngles = find(identities{m} == person);
    % Find the picture of this person under the same light angle
    Yextra{1} = [Yextra{1} ; Ytr{m}(allAngles(indsToSelect(i)),:)];
end
%{
for i=1:size(Yextra{1},1)
    for d=1:numberOfDatasets
        subplot(1,numberOfDatasets,d)
        imagesc(reshape(Yextra{d}(i,:),height, width)), title(num2str(identities{d}(i))),colormap('gray');
    end
    pause
end
%}
