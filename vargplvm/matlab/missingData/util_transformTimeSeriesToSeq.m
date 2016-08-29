% UTIL_TRANSFORMTIMESERIESTOSEQ: Given a timeseries, transform it in an
% auto-regressive dataset.
% DESC Given multivariate timeseries Y = [y1, y2, ...) where yi is a
% vector, create a dataset of input output pairs of the form:
% ([y1, y2, ..., y_L], y_{L+1})
% ([y2, y3, ..., y_{L+1}], y_{L+2})
% ...
% where L is the given timeWindow.
% 
% All the inputs are stored in X and corresponding outputs in Ynew.
%
% SEEALSO: util_transformSeqToTimeSeries.m, demKstepAhead.m
% 
function [X, Ynew] = util_transformTimeSeriesToSeq(Y, timeWindow)

Ntr = size(Y,1);
D = size(Y, 2);

blocksNumber = Ntr - timeWindow;
X = zeros(blocksNumber, timeWindow*D);
Ynew = zeros(blocksNumber, D);
for i=1:blocksNumber
    tmp =Y(i:i+timeWindow-1,:)';
    X(i,:) = tmp(:)';
    Ynew(i,:) = Y(i+timeWindow,:);
end