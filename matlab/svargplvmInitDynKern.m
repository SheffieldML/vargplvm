% SVARGPLVMINITDYNKERN Initialize a dynamics kernel for shared var.gplvm
% DESC Initialize a dynamics kernel for a shared variational gplvm model.
% VARGPLVM


% ATTENTION: For the gradients we assume that the base kernel (rbf,
% matern etc) must be the FIRST one and if a second base kernel
% (not white or bias) exist must be the LAST one!!!!!!!!!!!!!!
if isfield(kern,'comp')
    fprintf('# Dynamics Kernel initialization: \n')
    kernFound = 0;
    for k=1:numel(kern.comp)
        type = kern.comp{i}.type;
        if strcmp(type, 'rbfperiodic') || strcmp(type, 'rbfperiodic2')
            if exist('periodicPeriod')
                kern.comp{k}.period = periodicPeriod;
                kern.comp{k}.factor = 2*pi/periodicPeriod;
            end
            fprintf(1,'\t # periodic period: %d\n',kern.comp{k}.period);
        elseif strcmp(type, 'whitefixed')
            if ~exist('whiteVar')
                whiteVar = 1e-6;
            end
            kern.comp{k}.variance = whiteVar;
            fprintf(1,'\t # fixedwhite variance: %d\n',whiteVar);
        elseif strcmp(type, 'white')
            if ~exist('whiteVar')
                %     whiteVar = 1e-4; % Some models have been trained
                %     with this!!
                whiteVar = 0.1;
            end
            fprintf(1,'\t # white variance: %d\n',whiteVar);
            kern.comp{k}.variance = whiteVar; % Usual values: 1e-1, 1e-3
        elseif strcmp(type, 'bias')
            if exist('biasVar')
                kern.comp{k}.bias = biasVar;
                fprintf('\t # bias variance: %d \n', biasVar);
            end
        end
        % The following is related to the expected number of
        % zero-crossings.(larger inv.width numerator, rougher func)
        if strcmp(type,'rbfperiodic') || strcmp(type,'rbfperiodic2') || strcmp(type,'rbf') || strcmp(type,'matern32')
            kern.comp{k}.inverseWidth = optionsDyn{i}.inverseWidth./(((max(t{i})-min(t{i}))).^2);
            kern.comp{k}.variance = 1;
            % This is a bit hacky: if this is the second time an
            % rbf, or rbfperiodic or... kernel is found, then the
            % second variance can be initialised to be smaller
            if kernFound
                if ~exist('secondVarInit')
                    kern.comp{k}.variance = 0.006;
                else
                    kern.comp{k}.variance = secondVarInit;
                end
                fprintf('\t # Second variance initialized to %d \n',kern.comp{k}.variance);

            end
            kernFound = k;
        end
    end
end