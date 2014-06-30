% SVARGPLVM_INIT Initialise options. If the field name of 'defaults' already exists as a
% variable, the globalOpt will take this value, otherwise the default one.
%
% COPYRIGHT: Andreas C. Damianou, 2012, 2013
% VARGPLVM

if ~exist('globalOpt')
    defaults = util_optCreate({'optionsType', 'svargplvm'}); 
    
    %% The final options structure will contain default values, apart for whatever
    % is already set as a variable
    fnames = fieldnames(defaults);
    for i=1:length(fnames)
        if ~exist(fnames{i}, 'var')
            globalOpt.(fnames{i}) = defaults.(fnames{i});
        else
            globalOpt.(fnames{i}) = eval(fnames{i});
        end
    end
    clear('defaults', 'fnames');
    
end