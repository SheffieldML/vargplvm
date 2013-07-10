% VARGPLVMCREATETOYDATA Create a trivial dataset for testing var. GPLVM
% VARGPLVM

t = linspace(0,4*pi,100);
% The original signal will be a cosine, a sine and a squared cosine.
Z1 = cos(t)';
Z2 = sin(t)';
Z3= (cos(t)').^2;

% Store the original signals. Normally, the corresponding latent functions
% found for these signals will not be ordered, i.e. Z1 does not necessarily
% correspond to x_1(t). But with fixed random seeds it seems that the following
% correspondence is done Z1 ->x_3(t), Z2->x_2(t), X3 ->x_1(t)
if dynamicsUsed
    Z{1} = Z3; Z{2} = Z2; Z{3} = Z1;
else
    Z{1} = Z1; Z{2} = Z2; Z{3} = Z3;
end

% Scale and center data
bias_Z1 = mean(Z1);
Z1 = Z1 - repmat(bias_Z1,size(Z1,1),1);
scale_Z1 = max(max(abs(Z1)));
Z1 = Z1 ./scale_Z1;

bias_Z2 = mean(Z2);
Z2 = Z2 - repmat(bias_Z2,size(Z2,1),1);
scale_Z2 = max(max(abs(Z2)));
Z2 = Z2 ./ scale_Z2;

Z3 = Z3 - repmat(mean(Z3),size(Z3,1),1);%
Z3 = Z3 ./ (max(max(abs(Z3))));%

noiseLevel = 0.5; % Default: 0.1
% Map 1-D to 10-D and add some noise
Z2p = Z2*rand(1,10);
Z2p = Z2p + noiseLevel.*randn(size(Z2p));
Z1p = Z1*rand(1,10);
Z1p = Z1p + noiseLevel.*randn(size(Z1p));
Z3p = Z3*rand(1,10);%
Z3p = Z3p + noiseLevel.*randn(size(Z3p));%

% Form dataset by concatenating the signals
Y = [Z1p Z2p Z3p];
t = t';


%{
% Alternatively
t = linspace(0,4*pi,100);
% The original signal will be a cosine, a sine and a squared cosine.
Z{1} = cos(t)';
Z{2} = sin(t)';
Z{3}= (cos(t)').^2;
Z{4} = cos(t)' + 0.3;
Z{5} = cos(t)' + 0.6;
Z{6} = cos(t)' + 0.9;
Z{7} = sin(t)' + 0.3;
Z{8} = sin(t)' + 0.6;
Z{9} = sin(t)' + 0.9;
Z{10} = (cos(t)').^2 + 0.3;

noiseLevel = 0.4; % Default: 0.1
Y=[];
for i=1:10
    bias = mean(Z{i});
    Z{i} = Z{i} - repmat(bias,size(Z{i},1),1);
    scale = max(max(abs(Z{i})));
    Z{i} = Z{i} ./scale;
    Z{i} = Z{i} + noiseLevel.*randn(size(Z{i}));
    Y = [Y Z{i}];
end
t = t';
%}
