function [model, globalOpt] = svargplvmReOptimiseModel(model, initVardistIters, itNo, newExperimentNo, varargin)

% SVARGPLVMREOPTIMISEMODEL Optimise an already optimised svargplvm model for more iterations
%
% SEEALSO : svargplvmOptimiseModel
%
% COPYRIGHT : Andreas C. Damianou, 2011

% SVARGPLVM

globalOptOrig = model.globalOpt;
globalOpt = globalOptOrig;

if nargin < 4
    newExperimentNo = globalOpt.experimentNo;
end
globalOpt.experimentNo = newExperimentNo;

globalOpt.initVardistIters = initVardistIters;
globalOpt.itNo = itNo;
model.globalOpt = globalOpt;

model = svargplvmOptimiseModel(model, varargin{:});


globalOptOrig.initVardistIters = [globalOptOrig.initVardistIters initVardistIters];
globalOptOrig.itNo = [globalOptOrig.itNo itNo];

globalOpt = globalOptOrig;
model.globalOpt = globalOpt;
