function [X, model] = MRDEmbed(Yall, varargin)

% MRDEMBED Embed given data split in M views Yall{i}, i=1,...,M
% into a lower dimensional structure with MRD
% This function is just a wrapper for demSvargplvmGeneric.m
%
%
% [X, model] = MRDEmbed(Y, varargin)
%
% varargin = {'field1', value1, 'field2', value2, ...} where fields
% correspond to fields of the global options structure globalOpt the
% defaults of which can be found in "defaults_svargplvm.m".
%
% Alternatively:
% MRDEmbed(Yall, 'globalOpt', globalOpt)
% or
% MRDEmbed(Yall, 'options', options)
% or
% MRDEmbed(Yall, 'globalOpt', globalOpt, 'optionsDyn', optionsDyn)
% or
% MRDEmbed(Yall, 'options', options, 'optionsDyn', optionsDyn)
%
% COPYRIGHT: Andreas C. Damianou, 2014
% 
% SEEALSO: demSvargplvmGeneric.m
%
% VARGPLVM


% Default global options values
globalOptDef = util_optCreate({'optionsType', 'svargplvm'});
if nargin > 1 || ~isempty(varargin)
    if strcmp(varargin{1}, 'globalOpt') && isstruct(varargin{2})
       fprintf('# MRDEmbed: Structure globalOpt given!\n')
       globalOpt = varargin{2};
       givenOptions = true;
    elseif strcmp(varargin{1}, 'options') && isstruct(varargin{2})
        fprintf('# MRDEmbed: Structure options given!\n')
        options = varargin{2};
        givenOptions = true;
    else
        % Whatever different than default is given, overwrite it as an option.
        % varargin should be given as {'field1', value1, 'field2', value2, ...}
        globalOpt = util_optCreate(globalOptDef, varargin);
        givenOptions = false;
    end
    if givenOptions
        if nargin > 3 && strcmp(varargin{3}, 'optionsDyn') && isstruct(varargin{4})
            optionsDyn = varargin{4};
        end
    end
end

% The demo creates and optimises "model". The existing globalOpt structure
% will be used to initialise the model.
demSvargplvmGeneric;
X = model.vardist.means;