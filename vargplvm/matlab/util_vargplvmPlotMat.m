function vargplvmPlotMat(X, inds)
% UTIL_VARGPLVMPLOTMAT
% Plot in subplots each of the dimensions of a matrix (or cell matrix) X.
% If argument inds is provited, only the according dimensions will be
% plotted
% VARGPLVM

if iscell(X)
    Xnew = [];
    for i=1:length(X)
        Xnew = [Xnew X{i}];
    end
    X = Xnew;
    
end

if nargin == 2
    Q = length(inds);
else
    Q = size(X,2);
    inds = 1:Q;
end


dim1 = ceil(sqrt(Q));
dim2 = ceil(Q / dim1);

figure
for i=1:Q
    subplot(dim1, dim2, i)
    plot(X(:,inds(i)))
    title(['dim' num2str(inds(i))])
end
