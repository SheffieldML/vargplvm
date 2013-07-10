function retainedScales = vargplvmRetainedScales(model, thresh)

% VARGPLVMRETAINEDSCALES A small wrapper that shows which scales are switched off after Bayesian optimisation with Bayesian GP-LVM
% VARGPLVM

if nargin < 2
    thresh = 0.005;
end

if strcmp(model.kern.type, 'cmpnd') && ...
        (strcmp(model.kern.comp{1}.type, 'rbfard2') || strcmp(model.kern.comp{1}.type, 'linard2'))
    s1 = model.kern.comp{1}.inputScales;
else
    s1 = model.kern.inputScales;
end
% Normalise values between 0 and 1
s1 = s1 / max(s1);

%  thresh = max(model.kern.comp{1}.inputScales) * 0.001;

retainedScales = find(s1 > thresh);
