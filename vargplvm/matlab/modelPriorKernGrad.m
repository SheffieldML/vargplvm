function [muqOrig SqOrig] = modelPriorKernGrad(dynModel)
% MODELPRIORKERNGRAD
% see vargpTimeDynamicsPriorKernGrad. This function here is just a wrapper.
% VARGPLVM

fhandle = str2func([dynModel.type 'PriorKernGrad']);
[muqOrig SqOrig] = fhandle(dynModel);

% if isfield(model, 'paramGroups')
%  g = g*model.paramGroups;
% end
