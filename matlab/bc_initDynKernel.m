function kern = bc_initDynKernel(kern, globalOpt, optionsDyn, indexInComp)
% BC_INITDYNKERNEL Initialise a dynamics kernel
% VARGPLVM

X = optionsDyn.t;

fprintf('# Initialising kernel %s', kern.type);
if nargin > 3
    fprintf(' in place %d of the comp structure.\n', indexInComp);
else
    fprintf('\n');
end

% This is a recursive function. The recursion ending condition is for a
% kernel to not be a compound.
if isfield(kern, 'comp')
    for i=1:length(kern.comp)
        if ~isfield(kern.comp{i}, 'index') || isempty(kern.comp{i}.index)
            curInds = 1:kern.comp{i}.inputDimension;
        else
            curInds = kern.comp{i}.index;
        end
        curOptionsDyn = optionsDyn;
        curOptionsDyn.t = X(:, curInds);
        kern.comp{i} = bc_initDynKernel(kern.comp{i}, globalOpt, curOptionsDyn, i);
    end
else
    % This is not a compound kernel, so we can initialise this signel
    % element.
    kernelType = kern.type;
    if strcmp(kernelType, 'white')
        kern.variance = 1e-2; % Usual values: 1e-1, 1e-3
    end
    
    if strcmp(kernelType, 'whitefixed')
        kern.variance = globalOpt.fixedwhiteVar;
        fprintf(1,'# fixedwhite variance: %d\n',globalOpt.fixedwhiteVar);
    end
    

    % The following is related to the expected number of
    % zero-crossings.(larger inv.width numerator, rougher func)
    if strcmp(kernelType,'rbf') || strcmp(kernelType,'matern32') || strcmp(kernelType,'matern52')
        % NEW
        % dst = dist2(model.m, model.m);
        % lb = min(dst(:));
        % ub = max(dst(:));
        % inv_width = 2 / (ub+lb);
        %kern.comp{1}.inverseWidth = inv_width;
        
        kern.inverseWidth = optionsDyn.inverseWidth./(((max(max(optionsDyn.t))-min(min(optionsDyn.t)))).^2);

        if ~globalOpt.mappingInitialisation && (exist('indexInComp') && indexInComp == 1)
            kern.variance = 1;
        end
    end
end

