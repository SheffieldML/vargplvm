function [idcs, nSmpls] = class_samples( lbls, nActive )

% CLASS_SAMPLES
% COPYRIGHT: Patrick Sauer, 2012
% VARGPLVM

    if size(lbls,1) ~= 1
        lbls = lbls';
    end
    
    clids  = unique(lbls);
    nEx    = zeros(1, numel(clids));
    nSmpls = zeros(1, numel(clids));
    idcs   = cell( 1, numel(clids));
    
    count = 1;
    for c = clids 
        tmp = find(lbls==c);
        idcs{count} = tmp;
        nEx(count)  = numel(tmp);
        nSmpTmp =  ceil(nEx(count)/length(lbls)*nActive);
        if nSmpTmp > nEx(count)
            nSmpTmp = nEx(count);
        end
        nSmpls(count) = nSmpTmp;
       
        count = count + 1;
    end
    
    if sum(nEx) < nActive
       error('There must be at least as many training samples as inducing variables.');
    end
    
    nRes = nActive - sum(nSmpls);
    if nRes > 0
        while nRes > 0
            [nSmplsSort, ind] = sort(nSmpls,'ascend');
            for count = ind
                if nSmpls(count) < nEx(count)
                    nSmpls(count) = nSmpls(count)+1;
                    nRes = nRes-1;
                end

                if nRes == 0
                    break;
                end
            end
        end
    else
       while nRes < 0
          [mx, ind] = max(nSmpls);
          nSmpls(ind) = nSmpls(ind) - 1;
          nRes = nRes + 1;
       end
    end
end
