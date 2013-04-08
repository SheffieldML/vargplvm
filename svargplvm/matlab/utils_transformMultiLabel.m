function [Ynew, lblsNew] = utils_transformMultiLabel(Y, lbls)

% UTILS_TRANSFORMMULTILABEL Similar to transformLabels but this allows to handle multilabel datasets (not tested...)
% 
% SEEALSO: transformMultiLabel
%
% COPYRIGHT: Andreas C. Damianou, 2012
% SHEFFIELDML

Ynew = [];
lblsNew = [];
L = size(lbls,2);
for i=1:size(Y,1)
    indx = find(lbls(i,:));
    for j=1:length(indx)
        Ynew = [Ynew; Y(i,:)];
        newLbl = zeros(1,size(lbls,2));
        newLbl(indx(j)) = 1;
        lblsNew = [lblsNew; newLbl];
    end
end