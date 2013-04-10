function model = vargplvmDynamicsUpdateStats(model)
% VARGPLVMDYNAMICSUPDATESTATS wrapper function which according to the type
% of the model dynamics calls the appropriate function to perform the
% precomputations needed for the dynamics model
% DESC 
% ARG model: the model that contains the dynamics for which the statistics
% must be updated
% RETURN model : the updated model
%
% COPYRIGHT Andreas C. Damianou, Michalis Titsias, Neil Lawrence, 2010-2011
%
% VARGPLVM

fhandle = str2func([model.dynamics.type 'UpdateStats']);
model = fhandle(model);
