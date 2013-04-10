function ll = modelVarPriorBound(model)
% MODELVARPRIORBOUND Wrapper function for the various types of the
% variational GPLVM bound when dynamics is used.
% FORMAT
% DESC provides a wrapper function for the variational GP-LVM, when there
% are dynamics. It takes a model structure and according to its dynamics' type it
% calls the appropriate function to calculate the bound that corresponds
% only to the dynamics part (a KL term). The bound is a lower Jensens bound for
% variational approximation.
%
% ARG model : the model structure which contains the dynamics one, for
% which the corresponding part of the bound is to be computed.
% RETURN ll : the term of the variational bound which corresponds to the 
% dynamics structure found in the model.
% 
% SEEALSO : vargplvmLogLikelihood, vargpTimeDynamicsVarPriorBound
%
% COPYRIGHT : Michalis K. Titsias, 2010-2011
% COPYRIGHT : Neil D. Lawrence, 2010-2011
% COPYRIGHT : Andreas C. Damianou, 2010-2011


% VARGPLVM

fhandle = str2func([model.dynamics.type 'VarPriorBound']);
ll = fhandle(model.dynamics);