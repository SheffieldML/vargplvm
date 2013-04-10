function [bestThreshold, bestError, totalError] = svargplvmOptimiseThreshold(ZpredOrig, ZtrueOrig)

% SVARGPLVMOPTIMISETHRESHOLD Optimise the "binarization" threshold for classification
% DESC The model predicts labels in the continuous space, although labels are binary. If labels are between -1 and 1, 
% a straightforward idea is to tansform all predictions <0 to 1, and the rest to -1. However, due to biases, we could
% do even better by learning the best threshold on the training (or validation) sets.
%
% SEEALSO : demClassificationGeneralSvargplvm,  svargplvmPredictions3
%
% COPYRIGHT : Andreas C. Damianou, 2011
%
% VARGPLVM

L = size(ZtrueOrig,2);
ss=0;

bestThreshold = zeros(1,size(ZpredOrig,2));
bestError = ones(1, size(ZpredOrig,2)) * 1000;

if min(min(ZtrueOrig)) == -1
    start = -1;
else
    start = 0;
end

start = start + 0.1;
for d=1:size(ZpredOrig,2)
    for threshold = [0.1:0.0002:0.9]
        Ztrue = ZtrueOrig(:,d);
        Zpred = ZpredOrig(:,d);
        Zpred(Zpred >= threshold) = 1;
        Zpred(Zpred < threshold) = 0;
        hammingLoss = sum(sum(abs(Zpred - Ztrue))) / (size(Zpred,1) * size(Ztrue,2));
        if hammingLoss < bestError(d)
            bestThreshold(d) = threshold;
            bestError(d) = hammingLoss;
        end
    end
    ss = ss+bestError(d);
end

totalError = ss / L;
%ss / |L| gives the total hammingLoss
