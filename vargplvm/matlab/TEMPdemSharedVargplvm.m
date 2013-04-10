% DEM...
% DESC ...
%
% COPYRIGHT: Andreas C. Damianou, Michalis K. Titsias, 2010 - 2011
% SEEALSO: ...
% VARGPLVM


%close all; initLatent='ppca'; dataSetNames = {}; toyDataCreate ='vargplvm';initVardistIters = 400; itNo = [100 200 1500 500 500 400];  demSharedVargplvm1
%close all; initLatent='pca'; dataSetNames = {}; toyDataCreate ='fols';initVardistIters = 170; itNo = 300; dataToKeep=50;  demSharedVargplvm1
%clear ;close all; initLatent='pca3'; dataSetNames = {}; toyDataCreate='humanPose';initVardistIters = 160; itNo = [200 200 200 200 200 200 200 200]; dataToKeep=418;  demSharedVargplvm1

%clear ;close all; experimentNo=1;initLatent='pca3'; dataSetNames = {};toyDataCreate='humanPose';initVardistIters = 180; itNo = [500 200 200 200 200 200 200 200 200];  demSharedVargplvm1

clear timeStamps; % in case it's left from a previous experiment

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

% Define constants (in a manner that allows other scripts to parametrize
% this one).
if ~exist('experimentNo') ,  experimentNo = 404;      end
if ~exist('itNo')         ,  itNo = [500 500];              end     % Default: 2000
if ~exist('indPoints')    ,  indPoints = 80;          end     % Default: 49
if ~exist('latentDimPerModel')    ,  latentDimPerModel = 3;          end
% Set to 1 to use dynamics or to 0 to use the standard var-GPLVM
if ~exist('dynUsed')      ,  dynUsed = 0;             end
if ~exist('initVardistIters'), initVardistIters = 100;      end
if ~exist('dynamicKern')   ,     dynamicKern = {'rbf', 'white', 'bias'}; end
if ~exist('mappingKern')   ,  mappingKern = {'rbfard2', 'white'}; end
% if ~exist('mappingKern'),  mappingKern = {'rbfard2', 'bias', 'white'}; end
% Set to 1 to tie the inducing points with the latent vars. X
if ~exist('fixInd')        ,     fixInd = 0;    end
% 0.1 gives around 0.5 init.covars. 1.3 biases towards 0.
if ~exist('vardistCovarsMult'),  vardistCovarsMult=1.3;                  end
% Set to empty value {} to work with toy data
if ~exist('dataSetNames')    ,    dataSetNames = {'USecon','NYSE2Small'};    end
if ~exist('invWidthMultDyn'),    invWidthMultDyn = 100;                     end
if ~exist('invWidthMult'),       invWidthMult = 5;                     end
if ~exist('initX'),     initX ='ppca';   end % That's for the dynamics initialisation
if ~exist('initLatent'),     initLatent ='ppca';   end % That's for the whole latent space init.
if ~exist('dataToKeep'), dataToKeep = 100; end
if ~exist('toyDataCreate'), toyDataCreate = 'vargplvm'; end
if ~exist('doPredictions'), doPredictions = 0; end

%%
%-- Toy (random) data (samples from GPs)
if isempty(dataSetNames) || (isstr(dataSetNames) && strcmp(dataSetNames,'toy'))
    fprintf('# Creating toy data of type %s...\n', toyDataCreate);
    toyData = 1;
    %toyDataCreate = 'vargplvm'; % options: 'vargplvm', 'human','sampling', 'fols'
    switch toyDataCreate
        case 'humanPose'
            load humanPose;
            %ind2 = floor(1:2:size(Y,1));
            %%%__silhouettes are whole images
            if exist('imageSilhouette') && imageSilhouette
                load sil_images
            end
            % Remove the 'drunk' sequence
            Y1 = Y(1:100,:);
            Y2 = Y(337:end,:);
            Y = [Y1; Y2];
            clear('Y1','Y2');
            Z1 = Z(1:100,:);
            Z2 = Z(337:end,:);
            Z = [Z1; Z2];
            clear('Z1','Z2');
            % Subsample
            ind2 = floor(1:2:size(Y,1));
            Ytoy{1} = Y(ind2,:);
            Ytoy{2} = Z(ind2,:);
            dataSetNames={'silhouette', 'pose'};
            % X_init = Xp;
            %mappingKern = {'linard2', 'white'};
            %mappingKern = {'rbfard2', 'white'};
            latentDim = 5; % Anything > 2 and < 10
            %xyzankurAnim(Z_test, 3);
        case 'fols'
            %             addpath('../../ncca/matlab/');
            %             addpath('../../sgplvm/matlab/');
            %             addpath('../../fgplvm(svn)/matlab/');
            %             dem_sgplvm_fols; % This creates the dataset
            %             clear model
            %             clear options
            %alpha = linspace(0,2*pi,100);
            alpha = linspace(0,4*pi,100);
            Z1 = cos(alpha)';
            Z2 = sin(alpha)';
            Z3= (cos(alpha)').^2;
            % Z3 = 2*cos(2*alpha)' + 2*sin(2*alpha)' ; %
            
            
            % Scale and center data
            bias_Z1 = mean(Z1);
            Z1 = Z1 - repmat(bias_Z1,size(Z1,1),1);
            scale_Z1 = max(max(abs(Z1)));
            Z1 = Z1 ./scale_Z1;
            
            bias_Z2 = mean(Z2);
            Z2 = Z2 - repmat(bias_Z2,size(Z2,1),1);
            scale_Z2 = max(max(abs(Z2)));
            Z2 = Z2 ./ scale_Z2;
            
            Z3 = Z3 - repmat(mean(Z3),size(Z3,1),1);%
            Z3 = Z3 ./ (max(max(abs(Z3))));%
            
            noiseLevel = 0.3; % Default: 0.1
            % Map 1-D to 10-D and add some noise
            Z2p = Z2*rand(1,10);
            Z2p = Z2p + noiseLevel.*randn(size(Z2p));
            Z1p = Z1*rand(1,10);
            Z1p = Z1p + noiseLevel.*randn(size(Z1p));
            
            Z3p = Z3*rand(1,10);%
            Z3p = Z3p + noiseLevel.*randn(size(Z3p));%
            
            % pca(Z2p) % This shows that it is actually a 1-D dataset
            
            % Y = [Z1p Z2p];
            % We like the numer of latent dims to be 2+numSharedDims, ideally 3. With
            % vargplvm we set Q=6 and expect the other 3 or 4 to be switched off.
            % [U,V] = pca(Y,6);
            % Xp = Y*V;
            % pca(Xp)
            
            %---
            numSharedDims = 5;
            Z1p(:,1:numSharedDims) = Z3p(:,1:numSharedDims);
            Z2p(:,1:numSharedDims) = Z3p(:,1:numSharedDims);
            bar(pca([Z1p Z2p]))
            % return
            %---
            Ytoy{1} = Z1p;
            Ytoy{2} = Z2p;
            %Ytoy{3}= Z3p;%
            dataSetNames={'fols_cos', 'fols_sin'};
            % X_init = Xp;
            %mappingKern = {'linard2', 'white'};
            mappingKern = {'rbfard2', 'white'};
            latentDim = 6; % Anything > 2 and < 10
           % indPoints = 80;
        case 'vargplvm'
            % In the 'vargplvm' case, we generate X from a big Y (which
            % contains signal from rbf in the first half, and from a matern32
            % in the second half of the dims). Then, we set manually the
            % lengthscales for model_i to favor the i partition of X, and
            % generate a new Y_i from the posterior, hopefully associated with
            % the partition X_i. In the end (after optim), we would expect the lengthscales
            % found to be like the ones we manually set here.
            
            % Constants
            n = dataToKeep;
            D = 20;
            numSubModels = 2;
            latentDimPerModel = 3;
            numSharedDims = 0; % 2;
            %
            latentDim = latentDimPerModel * numSubModels + numSharedDims;
            x = linspace(-1, 1, n)';
            indPoints = round(2/3 * n);
            
            % Create toy Y as a combination of an rbf and a matern signal
            krn = 'rbf';
            invWidth = 10;
            kern = kernCreate(x,krn);
            kern.inverseWidth = invWidth;
            Kshared = kernCompute(kern, x);
            Y = zeros(n, D);
            for d=1:round(D/2)
                Y(:,d) = real(gsamp(zeros(1, size(x, 1)), Kshared, 1))';
            end
            
            krn = 'matern32';
            invWidth = 10;
            kern = kernCreate(x,krn);
            kern.inverseWidth = invWidth;
            Kshared = kernCompute(kern, x);
            for d=round(D/2+1):D
                Y(:,d) = real(gsamp(zeros(1, size(x, 1)), Kshared, 1))';
            end
            
            % Learn a model and X that fit Y
            options = vargplvmOptions('dtcvar');
            options.kern = 'rbfardjit';
            options.numActive = indPoints;
            options.optimiser = 'scg';
            d = size(Y, 2);
            model = vargplvmCreate(latentDim, d, Y, options);
            model = vargplvmParamInit(model, model.m, model.X);
            model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));
            % Optional____
            model.initVardist = 1;
            model = vargplvmOptimise(model, 1, 100);
            model.initVardist = 0;
            model = vargplvmOptimise(model, 1, 300);
            %_____
            % Create new Y_i by manipulating the ard parameters
            % latentIndices = randperm(size( model.vardist.means,2));
            latentIndices = 1:size(model.vardist.means,2);
            allIndices = 1:size(latentIndices,2);
            startVal = 1;
            endVal = numSharedDims;
            sharedIndices = latentIndices(startVal:endVal);
            dataSetNames = {};
            for i=1:numSubModels
                modelTemp = model;
                startVal = endVal+1;
                endVal = startVal + latentDimPerModel - 1;
                curIndices = [latentIndices(startVal:endVal) sharedIndices];
                switchOffIndices = setdiff(allIndices, curIndices);
                modelTemp.kern.inputScales(switchOffIndices) = 0;
                modelTemp.kern.comp{1}.inputScales(switchOffIndices) = 0;
                modelTemp.kern.comp{1}.inputScales(sharedIndices) = (1/2)* modelTemp.kern.comp{1}.inputScales(sharedIndices); % Optional
                modelTemp.kern.inputScales(sharedIndices) = (1/2)* modelTemp.kern.inputScales(sharedIndices); % Optional
                figure, bar(modelTemp.kern.inputScales)
                Ytoy{i} = vargplvmPosteriorMeanVar(modelTemp, modelTemp.vardist.means);
                dataSetNames = {dataSetNames{:}, ['datasetNo' num2str(i)]};
            end
            numberOfDatasets = length(dataSetNames);
            clear('model','modelTemp','options','x','D','n','Y','kern','Kshared','d');
            %{
        case 'sampling' % DOESN'T WORK
            %-- 2nd way: Generate (sample) X from N(0,I), then arbitrarily split X into
            % X_1, X_2,..., X_n, X_s  and fit a GP to then find the posteriors
            % p(Y_1 | X_1, Xs), ..., p(Y_n | X_n, X_s), ... , p(Y_s | X_s). Then, use
            % these Y's as datasets and see if you can get the X's back!!
            n = dataToKeep; % N
            indPoints = round((2/3)*n);
            subModelsNum = 2;
            numSharedDims = 2;
            Q = subModelsNum * (latentDimPerModel+numSharedDims);
            Xsamp = zeros(n,Q);
            for i=1:n
                Xsamp(i,:) = randn(1,Q);
            end
            latentIndices = randperm(size(Xsamp,2));
            startVal = 1;
            endVal =  numSharedDims;
            sharedIndices = latentIndices(startVal:endVal);
            for i=1:subModelsNum
                startVal = endVal+1;
                endVal = startVal + latentDimPerModel;
                X_i = [Xsamp(startVal:endVal) Xsamp(sharedIndices)];
                % Fit a gp
                % ...
                % Ytoy{i} =
            end
            %}
        case 'human'
            load('../../sgplvm/matlab/nccaDemoData.mat')
            Ytoy{1} = Y_train;
            Ytoy{2} = Z_train;
            numberOfDatasets=2;
        otherwise
            %-- For private signals
            dataSetNames = {'rand1','rand2','rand3','rand4'};
            kernels = {'rbf', 'rbf', 'ou','ou'};
            n = dataToKeep; % N
            indPoints = round((2/3)*n);
            dims = {40,40,40,40};
            invWidths = {10,10,15,15};
            numberOfDatasets = length(dataSetNames);
            x = linspace(-1, 1, n)';
            
            %-- For shared signal
            krn = 'rbf';
            invWidth = 100;
            kern = kernCreate(x,krn);
            kern.inverseWidth = invWidth;
            Kshared = kernCompute(kern, x);
            maxD = max(cell2mat(dims));
            Yshared = zeros(n, maxD);
            for d=1:maxD
                Yshared(:,d) = real(gsamp(zeros(1, size(x, 1)), Kshared, 1))';
            end
            
            %-- Build the datasets by adding the shared with the private signals
            for i=1:numberOfDatasets
                Ytoy{i} = zeros(n, dims{i});
                krn = kernels{i};
                for d=1:dims{i}
                    kern = kernCreate(x, krn);
                    kern.inverseWidth = invWidths{i};
                    K = kernCompute(kern,x);
                    Ytoy{i}(:,d) = real(gsamp(zeros(1, size(x, 1)), K, 1))';
                end
                % Add the shared signal
                Ytoy{i} = Ytoy{i} + Yshared(:,1:dims{i});
            end
            
            %-- We can add one extra dataset: the shared one
            Ytoy{i+1} = Yshared;
            numberOfDatasets = numberOfDatasets + 1;
            dataSetNames = {dataSetNames{:}, 'shared'};
            
            clear('n','K','x','kernels','krn','kern','d');
    end
    if ~exist('numberOfDatasets')
        numberOfDatasets = length(Ytoy);
    end
    %     for i=1:numberOfDatasets
    %         figure
    %         plot(Ytoy{i})
    %     end
else
    toyData = 0;
    numberOfDatasets = length(dataSetNames);
end
%--


%%

learnScales = 0;
mAll=[];

%-- Load datasets
for i=1:numberOfDatasets
    if toyData
        Y = Ytoy{i};
    else
        Y = lvmLoadData(dataSetNames{i});
        if dataToKeep ~= -1 && dataToKeep <= size(Y,1)
            Y = Y(1:dataToKeep,:);
        end
    end
    dims{i} = size(Y,2);
    N{i} = size(Y,1);
    indTr = 1:N{i};
    indTs = setdiff(size(Y,1), indTr);
    Ytr{i} = Y(indTr,:); %Yts = Y(indTs,:);
    t{i} = linspace(0, 2*pi, size(Y, 1)+1)'; t{i} = t{i}(1:end-1, 1);
    timeStampsTraining{i} = t{i}(indTr,1); %timeStampsTest = t(indTs,1);
    d{i} = size(Ytr{i}, 2);
end

%%%%%% TEMP: N{i}'s must be the same!!
for i=2:numberOfDatasets
    if N{i} ~= N{i-1}
        error('The number of observations in each dataset must be the same!');
    end
end
%%%%


%-- Options for the models
for i=1:numberOfDatasets
    % Set up models
    options{i} = vargplvmOptions('dtcvar');
    options{i}.kern = mappingKern; %{'rbfard2', 'bias', 'white'};
    %indPoints = 80; %%%%%
    options{i}.numActive = indPoints;
    options{i}.optimiser = 'scg';
    
    % !!!!! Be careful to use the same type of scaling and bias for all
    % models!!!
    
    % scale = std(Ytr);
    % scale(find(scale==0)) = 1;
    %options.scaleVal = mean(std(Ytr));
    options{i}.scaleVal = sqrt(var(Ytr{i}(:)));
    if fixInd
        options{i}.fixInducing=1;
        options{i}.fixIndices=1:size(Ytr{i},1);
    end
end


%-- Create the normalised version of the datasets and concatenate
%!!!! Doesn't work if Y{i}'s have different sizes!!
for i=1:numberOfDatasets
    % Compute m, the normalised version of Ytr (to be used for
    % initialisation of X)
    bias = mean(Ytr{i});
    scale = ones(1, d{i});
    
    if(isfield(options{i},'scale2var1'))
        if(options{i}.scale2var1)
            scale = std(Ytr{i});
            scale(find(scale==0)) = 1;
            if(learnScales)
                warning('Both learn scales and scale2var1 set for GP');
            end
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
     X_init{1} = ppcaEmbed(m{1},12);
     X_init{2} = ppcaEmbed(m{2},5);
     X_init = [X_init{1} X_init{2}];
elseif strcmp(initLatent,'ppcaConcatenate')
    initFunc = str2func([initX 'Embed']);
    X_init = initFunc(mAll, latentDimPerModel * numSubModels + numSharedDims);
elseif strcmp(initLatent, 'pca')
    % We like the numer of latent dims to be 2+numSharedDims, ideally 3. With
    % vargplvm we set Q=6 and expect the other 3 or 4 to be switched off.
    [U,V] = pca(mAll,latentDim);
    X_init = mAll*V;
elseif strcmp(initLatent, 'pca2')
    % We like the numer of latent dims to be 2+numSharedDims, ideally 3. With
    % vargplvm we set Q=6 and expect the other 3 or 4 to be switched off.
    [U,V] = pca(m{1},latentDim);
    X_init{1} = m{1}*V;
    [U,V] = pca(m{2},latentDim);
    X_init{2} = m{2}*V;
    X_init = [X_init{1} X_init{2}];
elseif strcmp(initLatent, 'pca3')
    clear mAll
    % We like the number of latent dims to be 2+numSharedDims, ideally 3. With
    % vargplvm we set Q=6 and expect the other 3 or 4 to be switched off.
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
latentDim = size(X_init,2);

% Free up some memory
clear('Y','m')



%-- Create the sub-models: Assume that for each dataset we have one model.
% This can be changed later, as long as we find a reasonable way to
% initialise the latent spaces.
for i=1:numberOfDatasets
    %---- Here put some code to assign X to the global common X which must
    % be created by doing pca in the concatenation of Y's...After this
    % point, model{i}.X will be the same for all i's. TODO...
    fprintf(1,'# Creating the model...\n');
    options{i}.initX = X_init;
    model{i} = vargplvmCreate(latentDim, d{i}, Ytr{i}, options{i});
    model{i}.X = X_init; %%%%%%%
    model{i} = vargplvmParamInit(model{i}, model{i}.m, model{i}.X);
    model{i}.X = X_init; %%%%%%%
    
    inpScales = invWidthMult./(((max(model{i}.X)-min(model{i}.X))).^2); % Default 5
    
    %inpScales(:) = max(inpScales); % Optional!!!!!
    
    model{i}.kern.comp{1}.inputScales = inpScales;
    
    if strcmp(model{i}.kern.type, 'rbfardjit')
        model{i}.kern.inputScales = model{i}.kern.comp{1}.inputScales;
    end
    params = vargplvmExtractParam(model{i});
    model{i} = vargplvmExpandParam(model{i}, params);
    
    
    model{i}.vardist.covars = 0.5*ones(size(model{i}.vardist.covars)) + 0.001*randn(size(model{i}.vardist.covars));
    
    %-------- Add dynamics to the model -----
    if dynUsed
        fprintf(1,'# Adding dynamics to the model...\n');
        optionsDyn{i}.type = 'vargpTime';
        optionsDyn{i}.t=timeStampsTraining{i};
        optionsDyn{i}.inverseWidth=invWidthMultDyn; % Default: 100
        optionsDyn{i}.initX = X_init; % initX; % That's for the dynamics
        
        kern = kernCreate(t{i}, dynamicKern); % Default: {'rbf','white','bias'}
        
        %-- Initialize each element of the compound kernel
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
        
        optionsDyn{i}.kern = kern;
        optionsDyn{i}.vardistCovars = vardistCovarsMult; % 0.23 gives true vardist.covars around 0.5 (DEFAULT: 0.23) for the ocean dataset
        
        % Fill in with default values whatever is not already set
        optionsDyn{i} = vargplvmOptionsDyn(optionsDyn{i});
        model{i} = vargplvmAddDynamics(model{i}, 'vargpTime', optionsDyn{i}, optionsDyn{i}.t, 0, 0,optionsDyn{i}.seq);
        
        fprintf(1,'# Further calibration of the initial values...\n');
        model{i} = vargplvmInitDynamics(model{i},optionsDyn{i});
        
        %  to also not learn the last kernel's variance
        if numel(kern.comp) > 1 && exist('learnSecondVariance') && ~learnSecondVariance
            fprintf(1,'# The variance for %s in the dynamics is not learned!\n',kern.comp{end}.type)
            model{i}.dynamics.learnSecondVariance = 0;
            model{i}.dynamics.kern.comp{end}.inverseWidth = model{i}.dynamics.kern.comp{1}.inverseWidth/10; %%% TEMP
        end
    end
    
    model{i}.beta=1/(0.01*var(model{i}.m(:)));
    prunedModelInit{i} = vargplvmPruneModel(model{i});
    %disp(model{i}.vardist.covars)
end


%modelInit = model;%%%TEMP

% TODO:
%--  Unify models into a structure
model = svargplvmModelCreate(model);
model.dataSetNames = dataSetNames;
model.initLatent = initLatent;
model.experimentNo = experimentNo;
%%---
capName = 'sharedVargplvm1';
capName(1) = upper(capName(1));
capName = []; %%%
modelType = model.type;
modelType(1) = upper(modelType(1));
fileToSave = ['dem' capName modelType num2str(experimentNo) '.mat'];
%%---


%%%
params = svargplvmExtractParam(model);
model = svargplvmExpandParam(model, params);
%    model = svargplvmOptimise(model, 1, 200);

%%


display = 1;
%%%% Optimisation
% do not learn beta and sigma_f for few iterations for intitialization
if initVardistIters ~=0
    model.initVardist = 1;
    model = svargplvmPropagateField(model,'initVardist', model.initVardist);
    fprintf(1,'# Intitiliazing the variational distribution...\n');
    model = svargplvmOptimise(model, display, initVardistIters); % Default: 20
    %fprintf(1,'1/b = %.4d\n',1/model.beta);
    
    modelInitVardist = model;
    model.initVardistIters=initVardistIters;
end

model.initVardist = 0;
model = svargplvmPropagateField(model,'initVardist', model.initVardist);

model.iters = 0;



% Optimise the model.
for i=1:length(itNo)
    iters = itNo(i); % default: 2000
    fprintf(1,'\n# Optimising the model for %d iterations (session %d)...\n',iters,i);
    model = svargplvmOptimise(model, display, iters);
    model.iters = model.iters + iters;
    % fprintf(1,'1/b = %.4d\n',1/model.beta);
    %fprintf(1,'# 1/b=%.4f\n var(m)=%.4f\n',1/model.beta, var(model.mOrig(:)));
    
    % Save model
    fprintf(1,'# Saving %s\n',fileToSave);
    %prunedModel = vargplvmPruneModel(model);
    % prunedModelTr = vargplvmPruneModel(modelTr);
    save(fileToSave, 'model', 'prunedModelInit');
end
%prunedModelTr = prunedModel;
%save(fileToSave, 'model', 'prunedModelInit', 'prunedModelTr');

% for i=1:numberOfDatasets
%     figure, bar(prunedModelInit{i}.kern.comp{1}.inputScales);
%     title(['Init scales for dataset ' num2str(i)]);
% end
for i=1:numberOfDatasets
    figure, bar(model.comp{i}.kern.comp{1}.inputScales);
    title(['Final scales for dataset ' num2str(i)]);
end

%%
%--------------------- PREDICTIONS --------------%
if ~doPredictions
    return
end
if model.numModels ~=2
    error('Predictions cannot be made (currently) with sharedVargplvm unless the number of submodels used is 2.');
end
%%

obsMod = 1; % one of the involved sub-models (the one for which we have the data)
infMod = setdiff(1:2, obsMod);

% Find the dimensions that are shared for obsMod and infMod
if ~exist('sharedDims')
    thresh = max(model.comp{obsMod}.kern.comp{1}.inputScales) * 0.001;
    retainedScales{obsMod} = find(model.comp{obsMod}.kern.comp{1}.inputScales > thresh);
    thresh = max(model.comp{infMod}.kern.comp{1}.inputScales) * 0.001;
    retainedScales{infMod} = find(model.comp{infMod}.kern.comp{1}.inputScales > thresh);
    sharedDims = intersect(retainedScales{obsMod}, retainedScales{infMod});
end

% Find X_* only for the shared dimensions (Xs*):
if ~exist('privateDims')
    privateDims{infMod} = setdiff(1:model.comp{obsMod}.q, sharedDims);
end

testOnTraining=0;

if testOnTraining
    i = size(model.comp{obsMod}.y,1); % last observation
    % Find p(X_* | Y_*): If Y* is not taken from the tr. data, then this step
    % must be an optimisation step of q(x*). X_* and Y_* here refer to the
    % spaces of the submodel obsMod.
    y_star = model.comp{obsMod}.y(i);
    x_star = model.comp{obsMod}.vardist.means(i,:);
    varx_star = model.comp{obsMod}.vardist.covars(i,:);
    [ind, distInd] = nn_class(model.X(:,sharedDims), x_star(:,sharedDims), 10, 'euclidean');
    
    ZpredMu = zeros(length(ind), size(Ytoy{infMod},2));
    ZpredSigma = zeros(length(ind), size(Ytoy{infMod},2));
    for k=1:length(ind)
        [ZpredMu(k,:), ZpredSigma(k,:)] = vargplvmPosteriorMeanVar(model.comp{infMod}, model.X(ind(k),:));
    end
else
    testInd = 25;
    Yts = Y_test(testInd,:); % Now this is a matrix
    for i=1:size(Yts,1)
        % initialize the latent points using the nearest neighbour
        % from the training data
        dst = dist2(Yts(i,:), model.comp{obsMod}.y);
        [mind, mini] = min(dst);
        
        Init(i,:) = model.vardist.means(mini,:);
        vardistx = vardistCreate(model.comp{obsMod}.vardist.means(mini,:), model.q, 'gaussian');
        vardistx.covars = model.comp{obsMod}.vardist.covars(mini,:);
        model.comp{obsMod}.vardistx = vardistx;
        display=1;
        iters = 100;
        [x_star, varx_star, modelUpdated] = vargplvmOptimisePoint(model.comp{obsMod}, vardistx, Yts(i,:), display, iters);%%%
        numberOfNN = 10;
        [ind, distInd] = nn_class(model.X(:,sharedDims), x_star(:,sharedDims), numberOfNN, 'euclidean');
        
        x = model.X(ind,:);
        mu = vargplvmPosteriorMeanVar(model.comp{infMod}, x);
        
        ZpredK{i} = mu;
        
        ZpredMu(i,:) = mu(1);
        
        %ZpredSigma(i,:) = sigma;
    end
     
    errsumFull = sum((ZpredMu - Z_test(testInd,:)).^2);
    errorFull = mean(errsumFull);
    
    if toyData && strcmp(toyDataCreate,'humanPose')
        for j = 1:1:length(testInd)
            handle = xyzankurVisualise(ZpredK{j}(1,:),1);
            for k = 2:length(ind)
              %  xyzankurModify(handle,ZpredK{j}(k,:));            handle = xyzankurVisualise(ZpredK{j}(1,:),1);
                handle = xyzankurVisualise(ZpredK{j}(k,:),k);
                title(['Pose'])
                pause;
            end
        end
    end
end

