function ret = svargplvmScales(method, model, varargin)

% SHEFFIELDMLSCALES A small utility to get or set the ARD weights (scales) for an svargplvm model.
%
% COPYRIGHT : Andreas C. Damianou, 2011

% SHEFFIELDML

if ~isfield(model, 'numModels') && isfield(model, 'M')
    model.numModels = model.M;
end

if strcmp(method, 'get')
    for i=1:model.numModels
        if strcmp(model.comp{i}.kern.type, 'rbfardjit')
            scales{i} = model.comp{i}.kern.inputScales;
        else
            scales{i} = model.comp{i}.kern.comp{1}.inputScales;
        end
    end
    ret = scales;
elseif strcmp(method, 'set')
    if nargin < 2
        error('Function requires at least 3 arguments with the set method')
    end
    scales = varargin{1};
    if length(scales) == 1
        scales = repmat(median(scales{1}), 1, model.q);
    end

    for i=1:model.numModels
        if strcmp(model.comp{i}.kern.type, 'rbfardjit')
            model.comp{i}.kern.inputScales = scales;
        else
            model.comp{i}.kern.comp{1}.inputScales = scales;
        end
    end
    ret = model;
end