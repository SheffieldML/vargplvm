function X_init = svargplvmInitLatentSpace(Ytr, d, options, initLatent, varargin)
% SVARGPLVMINITLATENTSPACE Initialise the latent space for a SVARGPLVM
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

latentDim = varargin{1};
latentDimPerModel = varargin{2};
numSharedDims = varargin{3};

mAll=[];
%-- Create the normalised version of the datasets and concatenate
%!!!! Doesn't work if Y{i}'s have different sizes!!
for i=1:length(Ytr)
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

% Clear some variables
clear('Y','Ytoy','bias','scale','ind2');


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
elseif strcmp(initLatent,'ppca')
    %-- Learn initialisation through PCA: Perform mappings from Y_i to X_i
    % and concatenate X_i's to augment the X's dimensionality as more
    % datasets are added.
    %     X_init = [];
    %     for i=1:numberOfDatasets
    %         initFunc = str2func([initX 'Embed']);
    %         X_init_cur = initFunc(m{i}, latentDimPerModel);
    %         X_init = [X_init X_init_cur];
    %     end
    
    %X_init{1} = ppcaEmbed(m{1},7);
    %X_init{2} = ppcaEmbed(m{2},3);
    
    X_init{1} = ppcaEmbed(m{1},latentDimPerModel);
    X_init{2} = ppcaEmbed(m{2},latentDimPerModel);
    
    
    X_init = [X_init{1} X_init{2}];
elseif strcmp(initLatent,'ppcaConcatenate')
    initFunc = str2func([initX 'Embed']);
    numSubModels = length(Ytr);
    X_init = initFunc(mAll, latentDimPerModel * numSubModels + numSharedDims);
elseif strcmp(initLatent, 'pca')
    % We like the numer of latent dims to be 2+numSharedDims, ideally 3. With
    % vargplvm we set Q=6 and expect the other 3 or 4 to be switched off.
    [U,V] = pca(mAll,latentDim);
    X_init = mAll*V;
elseif strcmp(initLatent, 'pca2')
    % Only for two submodels.
    [U,V] = pca(m{1},latentDim);
    X_init{1} = m{1}*V;
    [U,V] = pca(m{2},latentDim);
    X_init{2} = m{2}*V;
    X_init = [X_init{1} X_init{2}];
elseif strcmp(initLatent, 'pca3')
    clear mAll
    % Like pca initialisation but favour the first model compared to the
    % second (only for two submodels)
    try
        [U,V] = pca(m{1},7);
        X_init{1} = m{1}*V;
        m{1}=[];%%%
        [U,V] = pca(m{2},3);
        X_init{2} = m{2}*V;
        X_init = [X_init{1} X_init{2}];
    catch e
        if strcmp(e.identifier, 'MATLAB:nomem')
            fprintf('# !!! Warning: Not enough memory to initialise with PCA! Initialising with %s instead...\n',initX);
        end
        initFunc = str2func([initX 'Embed']);
        X_init{1} = initFunc(m{1}, 7);
        X_init{2} = initFunc(m{2},3);
        X_init = [X_init{1} X_init{2}];
    end
end
clear mAll
