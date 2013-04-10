function [cmat] = confMatrix( x, xref, labels_tr, labels_tst, nn )

% CONFMATRIX Create a confusion matrix for classification tasks.
%
% VARGPLVM

    totalLabels = length(unique(labels_tr));
    cmat = zeros(totalLabels,totalLabels);
    for i = 1:size(x,1)    
        dst = dist2(x(i,:), xref);
        [dsts, s] = sort(dst, 'ascend');
        sl = labels_tr(s(1:nn));
        % predicted class
        if nn > 1
            C =  mode(double(sl));
        else
            C = sl;
        end
        l = labels_tst(i);
        cmat(l, C) = cmat(l, C) + 1;
    end
end