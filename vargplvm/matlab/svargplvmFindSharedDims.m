function [sharedDims, privateDims] = svargplvmFindSharedDims(model, thresh, printOut, modalities)
% SVARGPLVMFINDSHAREDDIMS Find automatically the shared/private dimensions
% of the segmented latent space, based on some threshold
% DESC  Find automatically the shared/private dimensions of the segmented
% latent space, based on some given threshold for the value of the ARD
% weights
% 
% EXAMPLE: [sharedDims, privateDims] = svargplvmFindSharedDims(model,[],[],{[1 3] 2});
% 
% SEEALSO : svargplvmOptimiseModel, svargplvmShowScales
%
% COPYRIGHT: Andreas C. Damianou, 2012, 2013
%
% VARGPLVM

if nargin < 3 || isempty(printOut), printOut = false; end
if nargin < 2 || isempty(thresh), thresh = 0.005; end

if nargin < 4
    if length(model.comp) ~= 2
        error('method not implemented yet for > 2 submodels!')
    end
    obsMod = 1; % one of the involved sub-models (the one for which we have the data)
    infMod = setdiff(1:2, obsMod);
else
    obsMod = modalities{1};
    infMod = modalities{2};
end

for i=1:length(obsMod)
    if isfield(model.comp{obsMod(i)}.kern, 'comp')
        sObs{i} = model.comp{obsMod(i)}.kern.comp{1}.inputScales;
    else
        sObs{i} = model.comp{obsMod(i)}.kern.inputScales;
    end
end

if isfield(model.comp{infMod}.kern, 'comp')
    sInf = model.comp{infMod}.kern.comp{1}.inputScales;
else
    sInf = model.comp{infMod}.kern.inputScales;
end

% Normalise values between 0 and 1
for i=1:length(sObs)
    sObs{i} = sObs{i} / max(sObs{i});
end
sInf = sInf / max(sInf);

%  thresh = max(model.comp{obsMod}.kern.comp{1}.inputScales) * 0.001;

for i=1:length(sObs)
    retainedScales{obsMod(i)} = find(sObs{i} > thresh);
end
%thresh = max(model.comp{infMod}.kern.comp{1}.inputScales) * 0.001;
retainedScales{infMod} = find(sInf  > thresh);

sharedDims = intersect(retainedScales{obsMod(1)}, retainedScales{infMod});
for i=2:length(obsMod)
    sharedDims = intersect(sharedDims, retainedScales{obsMod(i)});
end

for i=1:length(retainedScales)
    privateDims{i} = setdiff(retainedScales{i}, sharedDims);
end

if printOut
    fprintf('# Shared dimensions:          [%s]\n', num2str(sharedDims))
    for i=1:length(retainedScales)
        fprintf('# Private dimensions model %d: [%s]\n', i, num2str(privateDims{i}))
    end
end