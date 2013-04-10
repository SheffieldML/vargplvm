function [labelErrors, labelSharingErrors, wrongIndices] = findLabelErrors(realLabels, predLabels)

% FINDLABELERRORS Find the number of errors in the predicted vs the real
% labels
% COPYRIGHT: Andreas C. Damianou, 2012
% SVARGPLVM

%%% WARNING: This function only works if there are no multi-label examples!

% First, transform any encoding to binary
% (e.g. if it is 1 -1 -1 -> 1 0 0)
lb = min(min(realLabels));
ub = max(max(realLabels));
realLabels(realLabels == lb) = 0;
realLabels(realLabels == ub) = 1;

predLabels(predLabels >= 0.5) = 1;
predLabels(predLabels < 0.5) = 0;

% %--- Old code
%     realUb = find(realLabels == ub);
%     predUb = predLabels(realUb);
%     % The number of times that a true label wasn't found
%     labelErrors = length(realUb) - sum(predUb); %equivalent to: length(find(~predUb))
%     
 % The number of times more than one label was found
 
 labelSharingErrors = sum((sum(predLabels,2)>1));
 
 %--
 realLabels2 = transformLabels(realLabels);
 predLabels2 = transformLabels(predLabels);
 uniqueLabels = unique(realLabels2);
 labelErrors = 0;
 for j=1:length(uniqueLabels)
     trueClassInd = find(realLabels2 == uniqueLabels(j));
     wrongIndices{j} = find(~(predLabels2(trueClassInd) == uniqueLabels(j)));
     labelErrors = labelErrors + length(wrongIndices{j});
 end

