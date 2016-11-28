function Yall = util_createMRD_toy_data()

alpha = linspace(0,4*pi,100);
% Scale and center data
Z1 = scaleData(cos(alpha)', 2);
Z2 = scaleData(sin(alpha)', 2);
Z3 = scaleData((cos(alpha)').^2, 2); % OR: 2*cos(2*alpha)' + 2*sin(2*alpha)'
noiseLevel = 0.1; % Default: 0.1
% Map 1-D to 10-D and add some noise
Z2p = Z2*rand(1,10);
Z2p = Z2p + noiseLevel.*randn(size(Z2p));
Z1p = Z1*rand(1,10);
Z1p = Z1p + noiseLevel.*randn(size(Z1p));
Z3p = Z3*rand(1,10);%
Z3p = Z3p + noiseLevel.*randn(size(Z3p));%
%---
numSharedDims = 5;
Z1p(:,1:numSharedDims) = Z3p(:,1:numSharedDims);
Z2p(:,1:numSharedDims) = Z3p(:,1:numSharedDims);
bar(pca([Z1p Z2p]))
title(sprintf('PCA scales in the final %d-dimensional dataset.', size(Z1p,2)+size(Z2p,2)))
%---
Yall{1} = Z1p;
Yall{2} = Z2p;
