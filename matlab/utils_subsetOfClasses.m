function [classInds, Y, labels, lbls, timeStamps] = utils_subsetOfClasses(Y, dataNo, labels, lbls, timeStamps)

% UTILS_SUBSETOFCLASSES  Out of a large dataset, only select a few datapoints (and corresponding labels), with constraints on the number
% of the total number of data or on the amount of data per class.
% SEEALSO: demClassificationGeneral
% COPYRIGHT: Andreas C. Damianou, 2012

% SVARGPLVM

if nargin<5
    timeStamps = [];
end
if nargin < 4
    lbls = [];
end
dataRest = 0;

if nargin > 2 
    uniqueLabels = unique(labels);
    numClasses = length(uniqueLabels);
    
    % We can either select a specific number of points per class, or a
    % total number of points. In the second case, the classes will be
    % almost but not exactly equally balanced, if mod(|Y|,dataPerClass ~=0)
    if iscell(dataNo) && length(dataNo) == 2
        if strcmp(dataNo{1}, 'dataPerClass')
            dataPerClass = dataNo{2};
        elseif strcmp(dataNo{1}, 'dataToKeep')
            dataToKeep = dataNo{2};
            dataPerClass = floor(dataToKeep / numClasses);
            dataRest = mod(dataToKeep, numClasses);
        end 
    else
        dataPerClass = dataNo;
    end
else
    error('Function requires at least three arguments.')
end

if dataPerClass ~= -1

    classInds = [];
    for c=1:numClasses
        curInds = find(labels == uniqueLabels(c));
        indPerm = randperm(length(curInds));
        % This is if we want a fixed but uneven number of elements for
        % class 1,2,...numClasses
        if length(dataPerClass) == numClasses
            curInds = curInds(indPerm(1:dataPerClass(c)));
        else
            curInds = curInds(indPerm(1:dataPerClass));
        end
        classInds = [classInds curInds];
    end
    if dataRest
        indTemp = setdiff(1:size(Y,1),classInds);
        indPerm = randperm(length(indTemp));
        indTemp = indTemp(indPerm(1:dataRest));
        classInds = [classInds indTemp];
    end
    
    Y = Y(classInds,:);
    if ~isempty(lbls)
        lbls = lbls(classInds,:);
    end
    labels = labels(classInds);
    if ~isempty(timeStamps)
        timeStamps = timeStamps(classInds);
    end
end