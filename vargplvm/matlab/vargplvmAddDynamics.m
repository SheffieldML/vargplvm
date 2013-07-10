function model = vargplvmAddDynamics(model, type, varargin)
% VARGPLVMADDDYNAMICS Add a dynamics structure to the model.
%
%	Description:
%
%	VARGPLVMADDDYNAMICS(MODEL, TYPE, ...) adds a dynamics structure to the
%	VARGPLVM.
%	 Arguments:
%	  MODEL - the model to add dynamics to.
%	  TYPE - the type of dynamics model to add in.
%	  ... - additional arguments to be passed on creation of the
%	   dynamics model.
%
% COPYRIGHT : Andreas C. Damianou, 2010-2011
% COPYRIGHT : Michalis K. Titsias, 2010-2011
% COPYRIGHT : Neil D. Lawrence, 2010-2011
%	
%
%	See also
%	MODELCREATE
% VARGPLVM

type = [type 'Dynamics'];
model.dynamics = modelCreate(type, model.q, model.q, model.X, varargin{:}); % e.g. vargpTimeDynamicsCreate
params = vargplvmExtractParam(model);
model.dynamics.nParams = length(modelExtractParam(model.dynamics));
model.nParams = model.nParams + model.dynamics.nParams;
model = vargplvmExpandParam(model, params);

% This will be moved to "initDynamics"
%med = median(model.vardist.covars(:));
%if med < 0.3 || med > 0.7
%   warning('!!! Median value of variational covariances is %1.2f.\n', med);
%end
