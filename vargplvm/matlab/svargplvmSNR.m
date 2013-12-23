function SNR = svargplvmSNR(model)

% SVARGPLVMSNR A small utility to display the Signal to Noise ratio of an optimised svargplvm model
% DESC The signal to noise ratio is defined for every modality, and it is shows the amount of  true signal
% vs the amount of "noise" which is assumed to have generated the data.
% SEEALSO : svargplvmOptimiseModel
%
% COPYRIGHT : Andreas C. Damianou, 2012

% VARGPLVM


for i=1:model.numModels
    if model.comp{i}.DgtN && isfield(model.comp{i}, 'mOrig') && ~isempty(model.comp{i}.mOrig)
        varData = var(model.comp{i}.mOrig(:));
    else
        varData = var(model.comp{i}.m(:));
    end
    SNR{i} = varData/(1/model.comp{i}.beta);
end