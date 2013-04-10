% FMRICREATEDATASET Create a dataset out of 4-D fMRI scans.
% DESC This script reads a 4-D NIFTI dataset of fMRI scans, processes it
% and stores it into a single matrix Y. The NIFTI dataset is assumed to
% exist in the path returned by localDatasetsDirectoryLarge with a default
% folder structure, i.e. appending:
% 'fmri/fmriDataFinal/fName/analyze/functional/functional4D.nii'.
% where fName is the folder name for each dataset. For classification, we
% also need the corresponding boldDelay3.txt file.
% The preprocessing done is the following:
% 1) Load the .nii file and serialise every 3-D scan in one row vector y_i corresponding to a specific time point t_i (and omit the very last image which is the mean).
% 2) Completely remove the datapoints (rows) of Y corresponding to the 'base' entries in the boldDelay3.txt (for ALL runs).
% 3) Calculate mean_rest and var_rest, i.e. the mean(Y(indices,:)) and var(Y(indices,:)) where 'indices' only takes into account the rows of Y corresponding to entries having the 'rest' label in the boldDelay3.txt (taking into account ALL runs).
% 4) Substitute every y_i with y_i = (y_i - mean_rest)/var_rest;
% 5) Completely remove the datapoints (rows) of Y corresponding to 'rest' entries of boldDelay3.txt
% 6) We are now left with groups of 4 rows, each group corresponding to a single "trial", no rest entries, no base entries. Replace each of these groups with the mean of the datapoints belonging to that.

% COPYRIGHT : Andreas C. Damianou, 2011
%
% SEEALSO : demFmriVargplvm1.m

% VARGPLVM

%---- Constants----
% maskPath='lc1ms_deskulled.img'
fName = '19740112LEZE_201011260900';
%fName = '19860611SUJN_201009211400';
delim = filesep;
maskPath = [localDatasetsDirectoryLarge 'fmri' delim 'fmriDataFinal' delim fName delim 'mask' delim 'lc1ms_deskulled.img'];


if ~exist('displ'), displ = 0; end
if ~exist('applyZscore'), applyZscore = 0; end
if ~exist('applyMask'), applyMask = 1; end
% Allows to further crop some dimensions which have 0 variance
if ~exist('thresholdOnVar'), thresholdOnVar = 0.8; end
if ~exist('discardRestPeriods'), discardRestPeriods = 0; end

% Number of datapoints to load (leave blank or empty to load all)
%if ~exist('N'), N = 1*460; end % All: 8*460

%if mod(N,460)
%    warning('The number of datapoints N must be an integer number of runs (i.e. k*460)');
%end


% nii=load_nii('functional4D.nii',1:N);
%{
if strcmp(computer, 'PCWIN')
    [a,p]=system('type NIFTI_pathWin.txt');
elseif strcmp(computer, 'GLNXA64') % server
    [a,p]=unix('cat NIFTI_path_server.txt');
else
    [a,p]=unix('cat NIFTI_path.txt');
end
%}
p = [localDatasetsDirectoryLarge 'fmri' delim 'fmriDataFinal' delim fName delim 'analyze' delim 'functional' delim 'functional4D.nii'];
disp('# Loading .nii file...')
if exist('N')
    nii=load_nii(p,1:N);
else
    nii = load_nii(p);
    % The '-1' is used to remove the last image which is the mean of all
    % the rest.
    N = size(nii.img,4)-1;
end
dimX = nii.original.hdr.dime.dim(2);
dimY = nii.original.hdr.dime.dim(3);
dimZ = nii.original.hdr.dime.dim(4);
dimAll = dimX * dimY * dimZ;
Y = zeros(N,dimAll);

%%
disp('# Serlializing data into a single matrix...')
% Store images in Y serialized, so that Y(1,:) = [nii.dimX' nii.dimY' nii.dimZ']
for i=1:N
    % Store current 3D image by columns
    curImg = nii.img(:,:,:,i);
    Y(i,:) = curImg(:)';
end
clear('curImg','nii');

if exist('maskPath') && applyMask
    disp('# Applying mask...')
    maskNii=load_nii(maskPath);
    if maskNii.hdr.dime.dim(1) ~= 3
        error('Mask is not a 3-D image!');
    end
    if dimX ~= maskNii.original.hdr.dime.dim(2) || dimY ~= maskNii.original.hdr.dime.dim(3) || dimZ ~= maskNii.original.hdr.dime.dim(4)
        error('Mask and data do not agree in dimensions!');
    end
    Ymask = maskNii.img(:)';
    
    %--
    indMask = find(Ymask);
    if thresholdOnVar
        thresh = 0; % mean(var(Y))/240;
        indVar = find(var(Y) > thresh);
        indNew = intersect(indMask, indVar);
        Ymask = zeros(1,size(Y,2));
        Ymask(indNew) = 1;
        Y = Y(:, indNew);
    else
        Y = Y(:,indMask);
    end
end

%% -- Now load the targets
disp('# Loading the targets...')
delims = find(p == delim);
pathTemp = p(1:delims(end-2));
pathTemp = [pathTemp 'behavioural' delim 'chunksTargets_boldDelay3.txt'];
[activity,runNumber]=textread(pathTemp','%s%s','delimiter',' ');
activity = activity(1:end-1);
runNumber = runNumber(1:end-1);
runNumber = str2num(cell2mat(runNumber));

% Remove 'base' entries:
disp('# Removing base entries...')
ind = find(strcmp(activity,'base'));
ind = setdiff(1:length(activity),ind);
Y = Y(ind,:);
activity = activity(ind);
runNumber = runNumber(ind);

ind = find(strcmp(activity,'rest'));
if applyZscore
    disp('# Applying Zscore based on rest periords...')
    %-- Z-score based on rest periods
    % Find the sample mean based on rest periods
    restMean = mean(Y(ind,:));
    restVar = var(Y(ind,:));
    % Slow approach but memory friendly
    for i=1:size(Y,1)
        Y(i,:) = (Y(i,:) - restMean)/restVar;
    end
end

% Now discard the 'rest' datapoints or summarize them in groups of 6:
restIndex = ind; clear('ind');
if ~discardRestPeriods
    Yrest = Y(restIndex, :);
    activityRest = activity(restIndex);
    runNumberRest = runNumber(restIndex);
end
indNotRest = setdiff(1:length(activity),restIndex);
Y = Y(indNotRest,:);
activity = activity(indNotRest);
runNumber = runNumber(indNotRest);
% Add the rest points in the end, after summarizing every 6 by taking the
% mean
if ~discardRestPeriods
    if mod(length(activityRest),6) ~= 0
        error('Rest periods are modelled as 6 consecutive time points...')
    end
    counter = 1;
    i=1;
    while counter<=size(Yrest,1)
        Yrest(i,:) = sum(Yrest(counter:counter+5,:))/6;
        activityRest(i) = activityRest(counter);
        runNumberRest(i) = runNumberRest(counter);
        i=i+1;
        counter = counter+6;
    end
    Yrest = Yrest(1:size(Yrest,1)/6,:);
    activityRest = activityRest(1:size(activityRest,1)/6,:);
    runNumberRest = runNumberRest(1:size(runNumberRest,1)/6,:);
end


% Finally, take the mean of each 4-datapoints:
% (slow approach but memory efficient):
disp('# Taking means of every 4 datapoints...')
counter = 1;
i=1;
while counter<=size(Y,1)
    Y(i,:) = sum(Y(counter:counter+3,:))/4;
    activity(i) = activity(counter);
    runNumber(i) = runNumber(counter);
    i=i+1;
    counter = counter+4;
end
Y = Y(1:size(Y,1)/4,:);
activity = activity(1:size(activity,1)/4,:);
runNumber = runNumber(1:size(runNumber,1)/4,:);

if ~discardRestPeriods
    Y = [Y; Yrest];
    activity = [activity; activityRest];
    runNumber = [runNumber; runNumberRest];
end

restIndex = find(strcmp(activity,'rest'));

% All 'animals':
animalIndex = find(strncmp(activity,'a-',2));

%%


if displ
    % Create an alternative dataset with Z direction fixed (so that you can see
    % the trajectories over time)
    fixedZdim = 15;
    YfixedZ = zeros(N, dimX*dimY);
    for i=1:N
        % Store current 3D image by columns
        curImg = nii.img(:,:,:,i);
        % For the alternative dataset
        curImg = nii.img(:,:,fixedZdim,i);
        YfixedZ(i,:) = curImg(:)';
    end
end


%-------

if displ
    % Plot image imNo for Z's dimensions dZ=1:end.
    % In the other subplot show the variance for ALL images for the
    % corresponding Z direction kept fixed
    imNo = 15;
    for dZ = 1:dimZ
        subplot(2,1,1)
        resStart = dimX * dimY * (dZ-1) +1;
        resEnd = resStart + dimX * dimY;
        imagesc(reshape(Y(imNo, resStart:resEnd-1), dimX, dimY))
        subplot(2,1,2)
        imagesc(reshape(var(Y(:,resStart:resEnd-1)),dimX,dimY))
        pause
    end
    
    close all
    dZ = fixedZdim;
    for i=1:N
        resStart = dimX * dimY * (dZ)+1;
        resEnd = resStart + dimX * dimY;
        imagesc(reshape(Y(i, resStart:resEnd-1), dimX, dimY))
        pause
    end
end
width = dimX;
height = dimY;
info.dataset = fName;
info.runs = length(unique(runNumber));
info.width = width; info.height = height; info.dimZ = dimZ;
info.date = date;
info.applyZscore = applyZscore;
info.applyMask = applyMask;
info.thresholdOnVar = thresholdOnVar;
info.discardRestPeriods = discardRestPeriods;

keep('Y','width','height','dimZ', 'activity','runNumber','animalIndex', 'info', 'restIndex');
