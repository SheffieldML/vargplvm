function lbls = transformLabels(labels)

% TRANSFORMLABELS A small utility which transform the label vector: if labels are in 1-of-K encoding it transforms it into a vector of
% integers and vice versa.
% SEEALSO: utils_transformMultiLabel
%
% COPYRIGHT: Andreas C. Damianou, 2012

% VARGPLVM

if size(labels,1) == 1 || size(labels,2) ==1
    N = max(size(labels,1), size(labels,2));    
    % The labels are not 1-K-encoding
    uniqueLabels = unique(labels);
    numClasses = length(uniqueLabels);
    lbls = zeros(N, numClasses);
    for c = 1:numClasses
        lbls(labels == uniqueLabels(c),c) = 1;
    end
else
    multilabel = false;
    N = size(labels,1);
    % This is in case we have e.g. 1 -1 -1 instead of 1 0 0 
    labels(labels==-1)=0; 
   % numClasses = size(labels,2);
    lbls = zeros(1,N);
    for n=1:N
        indx = find(labels(n,:));
        if length(indx) == 1
            lbls(1,n) = indx;
        else
            % WARNING!!!! When point n belongs to more than one classes,
            % this function will RANDOMLY select one!!!! This should be
            % fixed in future releases!
            permIndx = randperm(length(indx));
            lbls(1,n) = permIndx(1);
            multilabel = true;
        end
    end
    if multilabel
        warning('This is a multilabel dataset, transformLabels just kept one label per point randomly!')
    end
end

