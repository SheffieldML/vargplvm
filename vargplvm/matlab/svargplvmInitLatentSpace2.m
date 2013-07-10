function [X_init, m] = svargplvmInitLatentSpace2(Ytr, globalOpt, options)
% SVARGPLVMINITLATENTSPACE2 Initialise the latent space for a SVARGPLVM
% model
% FORMAT
% DESC Initialise the latent space for a SVARGPLVM model.
% ARG Ytr : A cell array containing the different datasets for the
% sub-models.
% ARG d : A cell array with the dimensionalities of the elements of Ytr
% ARG options : A cell array with each element being the options for the
% corresponding sub-model of the svargplvm model.
% ARG initLatent : How to initialise the latent space. Possible options include
% pca in the concatenated datasets, pca in each dataset and then
% concatenation, ncca etc.
% ARG varargin : Additional parameters, depending on the initialisation
% type.
% RETURN X_init : the initial latent points.
%
% SEEALSO : demSharedVargplvm1
%
% COPYRIGHT : Andreas C. Damianou, Carl Henrik Ek, 2011

% SVARGPLVM

latentDim = globalOpt.latentDim;
latentDimPerModel = globalOpt.latentDimPerModel;
%numSharedDims = globalOpt.numSharedDims;
initX = globalOpt.initX;
initLatent = globalOpt.initial_X;

mAll=[];
%-- Create the normalised version of the datasets and concatenate
%!!!! Doesn't work if Y{i}'s have different sizes!!
for i=1:length(Ytr)
    d{i} = size(Ytr{i},2);
    % Compute m, the normalised version of Ytr (to be used for
    % initialisation of X)
    bias = mean(Ytr{i});
    scale = ones(1, d{i});
    
    if(isfield(options{i},'scale2var1'))
        if(options{i}.scale2var1)
            scale = std(Ytr{i});
            scale(find(scale==0)) = 1;
            if(isfield(options{i}, 'scaleVal'))
                warning('Both scale2var1 and scaleVal set for GP');
            end
        end
    end
    if(isfield(options{i}, 'scaleVal'))
        scale = repmat(options{i}.scaleVal, 1, d{i});
    end
    
    % Remove bias and apply scale.
    m{i} = Ytr{i};
    for j = 1:d{i}
        m{i}(:, j) = m{i}(:, j) - bias(j);
        if scale(j)
            m{i}(:, j) = m{i}(:, j)/scale(j);
        end
    end
    
    mAll = [mAll m{i}]; % Concatenation (doesn't work if different sizes)
end


if size(globalOpt.initX,1) ~= 1 % Check if initial X is already given as a matrix
    X_init = globalOpt.initX;
else
    initFunc = str2func([initX 'Embed']);
    % %-- Create shared X:
    % initFunc = str2func([initX 'Embed']);
    % X = initFunc(mAll, latentDim);
    if ~isstr(initLatent)
        X_init = initLatent;
    elseif strcmp(initLatent, 'ncca')
        %-- Learn Initialisation through NCCA ( FOR TWO DATASETS only) %%%!!!
        if size(Ytr) ~= 2
            error('ncca initialization only when there are two datasets!');
        end
        [Xsy Xsz Xy Xz] = nccaEmbed(Ytr{1},Ytr{2},uint8([7 7]),uint8(1),uint8([2 2]),true);
        Xs = (1/2).*(Xsy+Xsz);
        X_init = [Xy Xs Xz]; % sizes: 2,1,2
        X_init = (X_init-repmat(mean(X_init),size(X_init,1),1))./repmat(std(X_init),size(X_init,1),1);
    elseif strcmp(initLatent,'separately')
        fprintf(['# Initialising the latent space with ' initX ' separately for each modality, with Q=['])
        if iscell(latentDimPerModel)
            for i=1:length(Ytr), fprintf('%d ', latentDimPerModel{i}); end
        else
            fprintf('%d', latentDimPerModel);
        end
        fprintf(']...\n')
        X_init = [];
        for i=1:length(Ytr)
            if iscell(latentDimPerModel)
                X_init_cur = initFunc(m{i},latentDimPerModel{i});
            elseif isscalar(latentDimPerModel)
                X_init_cur = initFunc(m{i},latentDimPerModel);
            else
                error('Unrecognised format for latentDimPerModel')
            end
            X_init = [X_init X_init_cur];
        end
    elseif strcmp(initLatent,'concatenated')
        fprintf(['# Initialising the latent space with ' initX ' after concatenating modalities in Q = %d ...\n'], latentDim)
        X_init = initFunc(mAll, latentDim);
    elseif strcmp(initLatent, 'custom')
        clear mAll
        % Like pca initialisation but favour the first model compared to the
        % second (only for two submodels)
        assert(length(Ytr)==2, 'Custom initialisation only for 2 submodels!')
        try
            fprintf(['# Initialising the latent space with ' initX ' separately, with (%d, %d) dims. for each modality...\n'],latentDimPerModel{1},latentDimPerModel{2})
            X_init{1} = initFunc(m{1},latentDimPerModel{1});
            if latentDimPerModel{2} ~= 0
                X_init{2} = initFunc(m{2},latentDimPerModel{2});
            else
                % For some aplications, we do not want embedding...(e.g. when
                % one of the modalities are the labels for classification)
                X_init{2} = m{2};
            end
            X_init = [X_init{1} X_init{2}];
        catch e
            if strcmp(e.identifier, 'MATLAB:nomem')
                warning(['Not enough memory to initialise with ' initX '! Initialising with PPCA instead...']);
            end
            initFunc = 'ppcaEmbed';
            X_init{1} = initFunc(m{1}, latentDimPerModel{1});
            X_init{2} = initFunc(m{2},latentDimPerModel{2});
            X_init = [X_init{1} X_init{2}];
        end
    else
        error('Unrecognised option for latent space initialisation.')
    end
end
clear mAll
