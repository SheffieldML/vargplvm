% SVARGPLVMPREPAREDATA Load and prepare different kinds of data for the shared variational GPLVM.
% DESC Load and prepare different kinds of data for the shared variational GPLVM.
%
% COPYRIGHT: Andreas C. Damianou, Carl Henrik Ek,  2011
% SEEALSO : demSharedVargplvm1
%
% VARGPLVM

%-- Toy (random) data (samples from GPs)
if isempty(dataSetNames) || (isstr(dataSetNames) && strcmp(dataSetNames,'toy'))
    fprintf('# Creating toy data of type %s...\n', toyDataCreate);
    toyData = 1;
    %toyDataCreate = 'vargplvm'; % options: 'vargplvm', 'human','sampling', 'fols'
    dataType = toyDataCreate;
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
            mappingKern = {'linard2', 'white'};
            % mappingKern = {'rbfard2', 'white'};
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
        case 'rbfOu'
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
        otherwise
            %...
    end
    if ~exist('numberOfDatasets')
        numberOfDatasets = length(Ytoy);
    end
    %     for i=1:numberOfDatasets
    %         figure
    %         plot(Ytoy{i})
    %     end
    Yall = Ytoy;
    clear('Ytoy');
else
    toyData = 0;
    if isstruct(dataSetNames)
        for i=1:length(dataSetNames)
            Y{i} = svargplvmLoadData(dataSetNames{i});
        end
        numberOfDatasets = length(dataSetNames);
    else
        [Yall, lbls] = svargplvmLoadData(dataSetNames);
        if ~isempty(lbls)
            height = lbls(1); width = lbls(2);
        end
        numberOfDatasets = length(Yall);
    end
end
%--

if isempty(dataType)
    if length(dataSetNames) ==1
        dataType = dataSetNames{1};
    end
end
