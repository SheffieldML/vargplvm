function Y = svargplvmMapModalitiesLow(balanceModalityDim, Y, originalDim, modalityMappingsAll)

if nargin < 4
    modalityMappingsAll = cell(1,length(Y));
end

if isscalar(balanceModalityDim) && ~iscell(balanceModalityDim)
    tmp = balanceModalityDim;
    balanceModalityDim = cell(1, length(Y));
    for i=1:length(Y)
        balanceModalityDim{i} = tmp;
    end
end

if sum(cell2mat(balanceModalityDim)) > 0
    maxDval = -Inf;
    % Find maximum dimensionality
    for i=1:length(Y)
        if size(Y{i},2) > maxDval
            maxDval = size(Y{i},2);
        end
    end
    
    % Do the mapping
    for i=1:length(Y)
        curD = originalDim(i); %size(Y{i},2);
        if curD < maxDval && balanceModalityDim{i}
            fprintf('# Mapping modality %d from %d to %d dimensions!\n', i, maxDval, curD)
            if isempty(modalityMappingsAll{i})
                % Fix seed for the random mapping so that it's reproducible
                curSeed = rng;
                rng(123);
                modalityMapping = rand(curD, maxDval);
                rng(curSeed);
                %--
            else
                modalityMapping = modalityMappingsAll{i};
            end
            Ytmp = zeros(size(Y{i},1), curD);
            reverseMapping = modalityMapping'*pdinv(modalityMapping*modalityMapping');
            for n = 1:size(Ytmp,1)
                Ytmp(n,:) = Y{i}(n,:)*reverseMapping;
            end
            Y{i} = Ytmp;
        end
    end
end

