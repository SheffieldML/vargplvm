function X = randomEmbed(Y, dims, numRange)

if nargin < 3
    numRange = [0,1];
end

if dims > size(Y,2)
    warning('Requested embedding in Q > D!')
end

X = (numRange(2)-numRange(1)).*rand(size(Y,1), dims) + numRange(1);
