% UTILS_OPTCREATE
% Create a structure of options given fieldName - value pairs in args as:
% args = {'field1', value1, 'field2', value,...};
% This function can be used as a template in different implementations.
% This can be done by:
% a) providing here an extra field {'optionsType', your_options_type}
% b) having an implementation of a function:
%  ['defaults_' your_options_type]
% which returns default values for this type of options
%
% COPYRIGHT: Andreas C. Damianou, 2012, 2013
% VARGPLVM
function options = util_optCreate(args1, args2)

options = [];

if nargin < 1
    return
end

if isstruct(args1)
    % In this case, the struct is initialised with another struct given in
    % arg1. Then, args2 (if present) can just alter specific fields or add.
    allfields = fieldnames(args1);
    for i=1:length(allfields)
        options.(allfields{i}) = args1.(allfields{i});
    end
    if nargin > 1
        for i = 1:2:length(args2)
            options.(args2{i}) = args2{i+1};
        end
    end
else 
    for i = 1:2:length(args1)
        options.(args1{i}) = args1{i+1};
    end
end



%--- This can change depending on the implementation needs ----------

%-- Set with default values whatever is not already set
% Default values
if isfield(options, 'optionsType')
    fhandle = str2func(['defaults_' options.optionsType]);
    defaults = fhandle();
else
    defaults = {};
end
%----------------------------------------------------------------------

for i = 1:2:length(defaults)
    if ~isfield(options, defaults{i})
        options.(defaults{i}) = defaults{i+1};
    end
end