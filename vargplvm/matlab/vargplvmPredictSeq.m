function [X, varX] = vargplvmPredictSeq(dynModel, t_star, seqInit)
% VARGPLVMPREDICTSEQ Predict the postions of a number of latent points.
% This function is not correct / complete and is not currently used.

% VARGPLVM


seq = dynModel.seq;
if ~exist('seqInit')
    seqInit = 1;
end
if seqInit ~= 1
    datFrom = seq(seqInit-1)+1;
else
    datFrom = 1;
end
datEnd=seq(seqInit);



N_star = size(t_star,1); % number of test points
K_ss = kernDiagCompute(dynModel.kern, t_star);
K_star = kernCompute(dynModel.kern, dynModel.t(datFrom:datEnd), t_star);

X = K_star' * dynModel.vardist.means(datFrom:datEnd,:); % mean

varX = zeros(N_star, dynModel.q); % initialize variances
for q=1:dynModel.q
    invLambda = 1./dynModel.vardist.covars((datFrom:datEnd),q); 
    Lq = chol(dynModel.Kt(datFrom:datEnd,datFrom:datEnd) + diag(invLambda))'; 
    vq = Lq \ K_star;
    varX(:,q) = K_ss - sum(vq .* vq)';
end
