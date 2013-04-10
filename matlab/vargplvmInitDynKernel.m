% VARGPLVMINITDYNKERNEL Initialise a compound dynamics kernel
% VARGPLVM

        if isfield(kern,'comp')
            fprintf('# Dynamics Kernel initialization: \n')
            kernFound = 0;
            for i=1:numel(kern.comp)
                type = kern.comp{i}.type;
                if strcmp(type, 'rbfperiodic') || strcmp(type, 'rbfperiodic2')
                    if exist('periodicPeriod')
                        kern.comp{i}.period = periodicPeriod;
                        kern.comp{i}.factor = 2*pi/periodicPeriod;
                    end
                    fprintf(1,'\t # periodic period: %d\n',kern.comp{i}.period);
                elseif strcmp(type, 'whitefixed')
                    if ~exist('whiteVar')
                        whiteVar = 1e-6;
                    end
                    kern.comp{i}.variance = whiteVar;
                    fprintf(1,'\t # fixedwhite variance: %d\n',whiteVar);
                elseif strcmp(type, 'white')
                    if ~exist('whiteVar')
                   %     whiteVar = 1e-4; % Some models have been trained
                   %     with this!!
                        whiteVar = 0.1;
                    end
                    fprintf(1,'\t # white variance: %d\n',whiteVar);
                    kern.comp{i}.variance = whiteVar; % Usual values: 1e-1, 1e-3
                elseif strcmp(type, 'bias')
                    if exist('biasVar')
                        kern.comp{i}.bias = biasVar;
                        fprintf('\t # bias variance: %d \n', biasVar);
                    end
                end
                % The following is related to the expected number of
                % zero-crossings.(larger inv.width numerator, rougher func)
                if strcmp(type,'rbfperiodic') || strcmp(type,'rbfperiodic2') || strcmp(type,'rbf') || strcmp(type,'matern32')
                    kern.comp{i}.inverseWidth = optionsDyn.inverseWidth./(((max(t)-min(t))).^2);
                    kern.comp{i}.variance = 1;
                    % This is a bit hacky: if this is the second time an
                    % rbf, or rbfperiodic or... kernel is found, then the
                    % second variance can be initialised to be smaller
                    if kernFound
                        if ~exist('secondVarInit')
                            kern.comp{i}.variance = 0.006;
                        else
                            kern.comp{i}.variance = secondVarInit;
                        end
                        fprintf('\t # Second variance initialized to %d \n',kern.comp{i}.variance);
                        
                    end
                    kernFound = i;
                end
            end
        end
