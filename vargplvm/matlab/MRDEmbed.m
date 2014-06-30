function [X, model] = MRDEmbed(Yall, varargin)

% MRDEMBED Embed given data split in M views Yall{i}, i=1,...,M
% into a lower dimensional structure with MRD
% This function is just a wrapper for demSvargplvmGeneric.m
%
% COPYRIGHT: Andreas C. Damianou, 2014
% 
% SEEALO: demSvargplvmGeneric.m
%
% VARGPLVM

% [X, model] = MRDEmbed(Y, varargin)
%
% varargin = {'field1', value1, 'field2', value2, ...} where fields
% correspond to fields of the global options structure globalOpt the
% defaults of which can be found in "defaults_svargplvm.m".


% Default global options values
globalOptDef = util_optCreate({'optionsType', 'svargplvm'});
if nargin > 1 || ~isempty(varargin)
    % Whatever different than default is given, overwrite it as an option.
    % varargin should be given as {'field1', value1, 'field2', value2, ...}
    globalOpt = util_optCreate(globalOptDef, varargin);
end

% The demo creates and optimises "model". The existing globalOpt structure
% will be used to initialise the model.
demSvargplvmGeneric;
X = model.vardist.means;