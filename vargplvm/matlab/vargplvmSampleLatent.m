% Sample layer automatically
function [X, mu, sigma] = vargplvmSampleLatent(model, dim, X, startingPoint)
if nargin <4 || isempty(startingPoint)
     % This point will be initially drawn. Then, we will sample and alter
     % only one if its dimensions.
    startingPoint = 1;
end

if nargin < 4 || isempty(X)
    Xorig = model.vardist.means;
    N = size(Xorig,1);
    xmin = min(Xorig(:,dim));
    xmax = max(Xorig(:,dim));
    df = xmax - xmin;
    xmin = xmin - 4*df/N; % also catch some points before xmin
    xmax = xmax + 4*df/N;
    x = linspace(xmin,xmax, 3*N); % this is the series of changes made in a specific dimension
    X = repmat(Xorig(startingPoint,:), length(x),1); % Just select some initial point
    X(:,dim) = x';
end

if nargout > 2
    [mu sigma] = vargplvmPosteriorMeanVar(model, X);
else
    mu = vargplvmPosteriorMeanVar(model, X);
end