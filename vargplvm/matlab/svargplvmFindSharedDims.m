function [sharedDims, privateDims] = svargplvmFindSharedDims(model, thresh, printOut)
% SVARGPLVMFINDSHAREDDIMS Find automatically the shared/private dimensions
% of the segmented latent space, based on some threshold
% DESC  Find automatically the shared/private dimensions of the segmented
% latent space, based on some given threshold for the value of the ARD
% weights
%
% SEEALSO : svargplvmOptimiseModel, svargplvmShowScales
%
% COPYRIGHT: Andreas C. Damianou, 2012, 2013
%
% VARGPLVM

if nargin < 3 || isempty(printOut), printOut = false; end
if nargin < 2 || isempty(thresh), thresh = 0.005; end

if length(model.comp) ~= 2
    error('method not implemented yet for > 2 submodels!')
end

obsMod = 1; % one of the involved sub-models (the one for which we have the data)
infMod = setdiff(1:2, obsMod);

if isfield(model.comp{obsMod}.kern, 'comp')
    s1 = model.comp{obsMod}.kern.comp{1}.inputScales;
else
    s1 = model.comp{obsMod}.kern.inputScales;
end
if isfield(model.comp{infMod}.kern, 'comp')
    s2 = model.comp{infMod}.kern.comp{1}.inputScales;
else
    s2 = model.comp{infMod}.kern.inputScales;
end

% Normalise values between 0 and 1
s1 = s1 / max(s1);
s2 = s2 / max(s2);

%  thresh = max(model.comp{obsMod}.kern.comp{1}.inputScales) * 0.001;

retainedScales{obsMod} = find(s1 > thresh);
%thresh = max(model.comp{infMod}.kern.comp{1}.inputScales) * 0.001;
retainedScales{infMod} = find(s2  > thresh);
sharedDims = intersect(retainedScales{obsMod}, retainedScales{infMod});

privateDims{1} = setdiff(retainedScales{1}, sharedDims);
privateDims{2} = setdiff(retainedScales{2}, sharedDims);

if printOut
    fprintf('# Shared dimensions:          [%s]\n', num2str(sharedDims))
    fprintf('# Private dimensions model 1: [%s]\n', num2str(privateDims{1}))
    fprintf('# Private dimensions model 2: [%s]\n', num2str(privateDims{2}))
end