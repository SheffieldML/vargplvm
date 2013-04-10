function model = svargplvmEqualizeScales(model, globalOpt)

% SVARGPLVMEQUALIZESCALES 
% COPYRIGHT: Andreas C. Damianou, 2012
% VARGPLVM

if nargin < 2
    globalOpt = model.globalOpt;
end

% Replace all scales with the median

scales = svargplvmScales('get',model);
conditionScales1=(max(scales{1}) - min(scales{1})/median(scales{1}) > 30);
conditionScales =  conditionScales1 && globalOpt.equalizeScales == 2;

if globalOpt.equalizeScales == 1
    fprintf('# Equalizing scales...\n')
    %newScales = repmat(median(scales{1}), 1, model.q);
    newScales = log(scales{1}) - 2*min(log(scales{1}));
    newScales(find(newScales < 1e-12)) = 1e-12;
    model = svargplvmScales('set', model, newScales);
elseif conditionScales
    fprintf('# Equalizing scales...\n')
    %newScales = repmat(median(scales{1}), 1, model.q);
    newScales = log(scales{1}) + 2*abs(min(log(scales{1})));
    newScales(find(newScales < 1e-5)) = 1e-5;
    model = svargplvmScales('set', model, newScales);
    model = svargplvmEqualizeScales(model, globalOpt); %%% Recursive!
elseif conditionScales1
    warning('Scales are very dissimilar!')
end