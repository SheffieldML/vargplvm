% UTIL_TRANSFORMSEQTOTIMESERIES The reverse function of
% UTIL_TRANSFORMTIMESERIESTOSEQ
%
% SEEALSO: util_transformTimeSeriesToSeq.m, demKstepAhead.m
function Ynew = util_transformSeqToTimeSeries(X, Y, timeWindow)

assert(size(X,1) == size(Y,1));

N = size(X,1) + timeWindow;
D = size(X,2) ./ timeWindow;

Ynew = zeros(N, D);

for i = 1:size(X,1)
    Ynew(i:i+timeWindow - 1, :) = reshape(X(i,:), D, timeWindow)';
end
Ynew(end,:) = Y(end,:);