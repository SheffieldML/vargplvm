function f = svargplvmObjective(params, model)

% SVARGPLVMOBJECTIVE Wrapper function for shared variational GP-LVM objective.
% FORMAT
% DESC provides a wrapper function for the shared variational GP-LVM, it
% takes the negative of the log likelihood, feeding the parameters
% correctly to the model.
% ARG params : the parameters of the variational GP-LVM model.
% ARG model : the model structure in which the parameters are to be
% placed.
% RETURN f : the negative of the log likelihood of the model.
% 
% SEEALSO : svargplvmModelCreate, svargplvmLogLikelihood, svargplvmExpandParam
%
% COPYRIGHT : Andreas C. Damianou, 2011

% VARGPLVM


model = svargplvmExpandParam(model, params);

f = -svargplvmLogLikelihood(model);
