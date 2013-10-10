function factor = transformedParamChainRule(transformType, transformIndex, params)
% TRANSFORMEDPARAMCHAINRULE Apply the appropriate multiplications coming
% from the chain rule for transformed parameters that are placed a prior
% DESC Apply the appropriate multiplications coming
% from the chain rule for transformed parameters that are placed a prior
% SEE ALSO: vargplvmParamPriorGradients
% COPYRIGHT: Andreas C. Damianou, 2013
% VARGPLVM

% assert(length(transformIndex) <= length(params));

switch transformType
    case 'exp'
        factor = params;
        notIndexed = setdiff(1:length(params), transformIndex);
        factor(notIndexed) = 1; % Do not alter the non-indexed parameters
    otherwise
        error('Not implemented for the selected transform')
end