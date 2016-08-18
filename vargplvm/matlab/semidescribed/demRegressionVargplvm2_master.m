% Fix seeds
if ~exist('rseed'), rseed = 1e4; end

randn('seed', rseed);
rand('seed', rseed);

if ~exist('dataSetName', 'var'), dataSetName = 'toy'; end
if ~exist('Ntr', 'var'), Ntr = 40; end % 60
if ~exist('D', 'var'), D = 12; end % 30
if ~exist('Q', 'var'), Q = 12; end % 10
if ~exist('QmissingPerc', 'var'), QmissingPerc = 0.2; end
if ~exist('Ninterp', 'var'), Ninterp = 40; end % 50
if ~exist('Nmissing', 'var'), Nmissing = 40; end % 30
if ~exist('outputNoise', 'var'), outputNoise = 0.001; end
if ~exist('numActive', 'var'), numActive = 25; end
if ~exist('iters', 'var'), iters = 500; end
if ~exist('initIters', 'var'), initIters = 1000; end
if ~exist('initSNR', 'var'), initSNR = 80; end
if ~exist('init_missing', 'var'), init_missing = 'zeros'; end
if ~exist('linearMapping', 'var'), linearMapping = false; end
if ~exist('resultsFile', 'var'), resultsFile = ''; end
if ~exist('displayIters', 'var'), displayIters = true; end
if ~exist('runBGPLVMmisFixed', 'var'), runBGPLVMmisFixed = false; end
if ~exist('runGP', 'var'), runGP = true; end
if ~exist('runBGPLVM', 'var'), runBGPLVM = true; end
if ~exist('runGPLVM','var'), runGPLVM = false; end
if ~exist('cmu_dim_output', 'var'), cmu_dim_output = 1:20; end
if ~exist('runBGPLVM_missing', 'var'), runBGPLVM_missing = true; end

N = Ntr + Ninterp +  Nmissing;

indExtrap = Ntr+Ninterp+1:N;
perm = randperm(Ntr + Ninterp + Nmissing);
indTr = sort(perm(1:Ntr));
indInterp = sort(perm(Ntr + 1:Ntr+Ninterp));
indMissing = sort(perm(Ntr+Ninterp+1:end));
if strcmp(dataSetName, 'toy')
    rep = 3;
    t = linspace(0, rep*2*pi, N)';
    
    kern1 = kernCreate(t, 'rbf');
    kern1.inverseWidth = 0.2;
    Kx = kernCompute(kern1, t);
    X = gsamp(zeros(N, 1), Kx, Q)';
    Xtr = X(indTr, :);
    [~, bias_Xtr, scale_Xtr] = scaleData(Xtr, 1);
    X = scaleData(X, 0, scale_Xtr, bias_Xtr);
    
    kern2 = kernCreate(X, 'rbf');
    kern2.inverseWidth = 0.5;
    Kf = kernCompute(kern2, X);
    Y = gsamp(zeros(N, 1), Kf, D)';
    
    Y = Y + sqrt(outputNoise).*randn(size(Y));
elseif strcmp(dataSetName, 'cmu35gplvm')
    %%% !!!! ADD YOUR PATH TO DATASETS (from GPmat) here!!!!!!!
    Y = lvmLoadData('cmu35gplvm', [],'/DATASETS0p1371/');
    % Leg inds are 8:14 and body is 21:50
    restInd = setdiff(1:size(Y,2), cmu_dim_output);
    X = Y(:, restInd);
    Y = Y(:, cmu_dim_output);
    Q = size(X,2);
    D = size(Y,2);
end
Ytr = Y(indTr, :);           Xtr = X(indTr,:);
Yinterp = Y(indInterp, :);   Xinterp = X(indInterp, :);
Ymissing = Y(indMissing, :); Xmissing = X(indMissing, :);

% Scale to variance 1 and zero mean, (inputs already scaled above for toy)
[Ytr, bias_Ytr, scale_Ytr] = scaleData(Ytr, 1);
Yinterp = scaleData(Yinterp, 0, scale_Ytr, bias_Ytr);
Ymissing = scaleData(Ymissing, 0, scale_Ytr, bias_Ytr);
if ~strcmp(dataSetName, 'toy')
    [Xtr, bias_Xtr, scale_Xtr] = scaleData(Xtr, 1);
    Xinterp = scaleData(Xinterp, 0, scale_Xtr, bias_Xtr);
    Xmissing = scaleData(Xmissing, 0, scale_Xtr, bias_Xtr);
end

% Now remove *dimensions* from X: ones indicate missing value
maskMissing = rand(size(Xmissing)) < QmissingPerc;
%sum(sum(maskMissing)) ./ prod(size(maskMissing))
XtrMissing = Xmissing;
XtrMissing(maskMissing) = NaN;
XtrMissingTmp = XtrMissing;
% A mask for the missing indices of the extended set with training and
% missing
maskMissingExt = [false(size(Xtr)); maskMissing];

%% Multiple linear regression with missing inputs

% Fit model and make prediction for test points
YPredLin = zeros(size(Yinterp));
b = cell(size(Yinterp,2),1);
YtrAll = [Ytr ; Ymissing];
for d=1:size(Yinterp,2)
    b{d} = regress(YtrAll(:,d), [Xtr; XtrMissing]);
    YPredLin(:,d) = Xinterp*b{d};
end

%{
% Plot results of fit to test points
t = 1:size(Yinterp,1);
for d=1:D
    plot(t, Yinterp(:,d), 'x--')
    hold on
    plot(t, YPredLin(:,d), 'ro-')
    pause
    hold off
end
%}

%% NN
miniTs=[];
for i=1:size(Yinterp,1)
    dst = dist2(Yinterp(i,:), Ytr);
    [~, miniTs(i)]=min(dst);
end
YPredNN = Ytr(miniTs,:);


YPredNNx = nan(size(YPredNN));
Xall = [Xtr; XtrMissing];
for i=1:size(Xinterp,1)
    minDist = Inf;
    minInd = -1;
    for j=1:size(Xall,1)
        inds = find(~isnan(Xall(j,:)));
        curDist = dist2(Xinterp(i,inds), Xall(j,inds));
        if curDist < minDist
            minInd = j;
            minDist = curDist;
        end
    end
    YPredNNx(i,:) = YtrAll(minInd,:);
end
%%
if runGP
    fprintf('\n\n# GP \n\n')
    optionsGP = gpOptions('dtcvar');
    optionsGP.numActive = numActive;
    if linearMapping
        optionsGP.kern = {'linard2','white','bias'};
    else
        optionsGP.kern =  {'rbfard2','bias','white'};
    end
    modelGP = gpCreate(Q, D, Xtr, Ytr, optionsGP);
    modelGP.beta = 1/((1/initSNR * var(modelGP.m(:))));
    modelGP.optimiseBeta = false;
    modelGP = gpOptimise(modelGP, displayIters, initIters);
    modelGP.optimiseBeta = true;
    modelGP = gpOptimise(modelGP, displayIters, sum(iters));
    YpredGP = gpPosteriorMeanVar(modelGP, Xinterp);
end
%% This is like the GP, since covars are fixed close to zero. But from here
%% we have access to the posterior q(X).
if runBGPLVM
    fprintf('\n\n\n# BGPLVM \n\n')
    
    options = vargplvmOptions('dtcvar');
    if linearMapping
        options.kern = {'linard2','white','bias'};
    else
        options.kern = 'rbfardjit';
    end
    options.numActive = numActive;
    options.optimiser = 'scg2';
    options.initSNR = initSNR;
    options.initX = Xtr;
    [~,~,~, model] = vargplvmEmbed(Ytr, Q, options, 0, 0);
    model.vardist.covars = 1e-10*model.vardist.covars;
    model.fixParamIndices = 1:2*model.N*model.q;
    model = vargplvmOptimiseModel(model, 0, 0, {ceil(initIters/2), 0}, displayIters);
    model = vargplvmOptimiseModel(model, 0, 0, {ceil(initIters/2), iters}, displayIters);
    YpredBGPLVM = modelPosteriorMeanVar(model, Xinterp);
end
%%
if runBGPLVM_missing
    fprintf('\n\n\n# BGPLVM with missing \n\n')
    
    options = vargplvmOptions('dtcvar');
    if linearMapping
        options.kern = {'linard2','white','bias'};
    else
        options.kern = 'rbfardjit';
    end
    options.numActive = numActive;
    options.optimiser = 'scg2';
    options.initSNR = initSNR;
    
    mini =[];
    for i=1:size(Ymissing,1)
        dst = dist2([Ymissing(i,:) XtrMissing(i, ~maskMissing(i,:))], [Ytr Xtr(:, ~maskMissing(i,:))]);
        [mind, mini(i)] = min(dst);
    end
    
    
    switch init_missing
        case 'NN'
            warning('That was commented out')
            XtrMissing(maskMissing) = Xtr(mini) + randn(Nmissing,1); %%%
        case 'posterior' %%% TODO
            % error('For this you need runBGPLVM flag active.')
            initQx.means = model.vardist.means(mini,:);
            initQx.covars = repmat(0.5, size(model.vardist.covars(mini,:))); %model.vardist.covars(mini,:);
            initQx.means(~maskMissing) = XtrMissing(~maskMissing);
            initQx.covars(~maskMissing) = 1e-6;
            [x_star_all, varx_star_all] = vargplvmPredictLatent(model, Ymissing, [], 0, 1000, 0, [], initQx);
            XtrMissing(maskMissing) = x_star_all(maskMissing);
        otherwise
            XtrMissing(maskMissing) = zeros(size(XtrMissing(maskMissing))) + randn(size(XtrMissing(maskMissing)))*0.00000001;
    end
    options.initX = [Xtr; XtrMissing];
    
    tries = 15;
    nanFlag = 1;
    while tries > 0 && nanFlag == 1
        [~,C]=kmeans_matlab(options.initX, options.numActive, 'EmptyAction', 'singleton','Start','uniform');
        nanFlag = isnan(sum(sum(C)));
        tries = tries - 1;
    end
    if nanFlag
        fprintf('!!! ERROR! Tried Kmeans 15 times but NaN was returned in all tries! \n') 
        return
    end
    
    %plot(options.initX, zeros(size(options.initX,1),1), 'x'); hold on; plot(C, zeros(size(C,1),1), 'Or')
    
    options.initX_u = C;
    
    [~,~,~, model_m] = vargplvmEmbed([Ytr ; Ymissing], Q, options, 0, 0);
    model_m.throwSNRError = false;
    model_m.vardist.covars(~maskMissingExt) = 1e-10;%*model_m.vardist.covars(~maskMissingExt);
    model_m.vardist.covars(maskMissingExt) = 1;
    if strcmp(init_missing, 'posterior')
        tmp = model_m.vardist.covars(Ntr+1:end,:);
        %tmp(maskMissing) = varx_star_all(maskMissing)*5;
        %--
        tmp2 = varx_star_all(maskMissing);
        tmp2(tmp2>1) = 0.99;
        tmp(maskMissing) = tmp2;
        %--
        model_m.vardist.covars(Ntr+1:end,:) = tmp;
    end
    model_m.kern.inputScales = ones(1, model_m.q)*(mean(model_m.kern.inputScales)) + randn(1, model_m.q)*0.01;
    model_m_first = model_m;
    
    %%
    % This model is if we train on fully observed data (like a GP), and
    % then infer the positions of new data and augment the model without
    % reoptimising it.
    model_BGPLVM_init = model_m;
    model_BGPLVM_init.kern.inputScales = model.kern.inputScales;
    model_BGPLVM_init.beta = model.beta;
    YpredBGPLVM_init = modelPosteriorMeanVar(model_BGPLVM_init, Xinterp);
    %%
    
    % ---
    %{
    figure
    xtmp = [Xtr; XtrMissingTmp];
    for i=1:Q
        plot(xtmp(:,i), 'ro--'); hold on;
        plot(model_m.vardist.means(:,i), 'bx');
        plot(model_m.vardist.means(:,i)+2*sqrt(model_m.vardist.covars(:,i)), 'g--')
        plot(model_m.vardist.means(:,i)-2*sqrt(model_m.vardist.covars(:,i)), 'g--')
        pause; hold off
    end
    %}
    % ---
    
    % We want to find the indices of the training/missing means and
    % covariances. Because the means and covariances are stored in the params
    % as means(:)' (i.e columnwise), we follow the next (inefficient) trick
    % and construct a matrix which repeats a specific number if the parameter
    % in that place is a training mean, missing ind. mean, training cov. or
    % missing ind cov. Here, missing ind means the indices of the XtrMissing
    % matrix, which doesn't have all elements missing.
    CC = [zeros(Ntr, Q); ones(Nmissing,Q)+maskMissing];
    DD = CC+10;
    cc = [CC(:)' DD(:)'];
    
    trMeansInd = find(cc==0);
    trCovInd = find(cc==10);
    miMeansIndPresent = find(cc==1);
    miMeansIndMissing = find(cc==2);
    miCovIndPresent = find(cc==11);
    miCovIndMissing = find(cc==12);
    
    miMeansInd = sort([miMeansIndPresent miMeansIndMissing]);
    miCovInd = sort([miCovIndPresent miCovIndMissing]);
    
    fixParamIndicesPresent = sort([trMeansInd trCovInd miMeansIndPresent miCovIndPresent]);
    % All var. means, all var. covars
    fixParamIndicesAll = sort([trMeansInd trCovInd miMeansInd miCovInd]);
    
    model_m.fixParamIndices = fixParamIndicesAll;
    
    % Initialising var. distr
    model_m_init = vargplvmOptimiseModel(model_m, 0, 0, {ceil(initIters/2), 0}, displayIters);
    model_m_init = vargplvmOptimiseModel(model_m_init, 0, 0, {ceil(initIters/2), 0}, displayIters);
    model_m_fixed = model_m_init;
    
    if runGPLVM, model_gplvm = model_m_fixed; end
    
    fprintf('\n\n# Release missing vardist params!\n\n');
    model_m_init.fixParamIndices = fixParamIndicesPresent; % Release missing vardist elements
    model_m_init = vargplvmOptimiseModel(model_m_init, 0, 0, {[ceil(initIters/4) ceil(initIters/4)], 0}, displayIters);
    model_m = vargplvmOptimiseModel(model_m_init, 0, 0, {0, iters}, displayIters);
    YpredBGPLVM_missing = modelPosteriorMeanVar(model_m, Xinterp);
    
    if runBGPLVMmisFixed
        fprintf('\n\n# Optimising the fixed model!\n\n');
        model_m_fixed = vargplvmOptimiseModel(model_m_fixed, 0, 0, {[ceil(initIters/4) ceil(initIters/4)], 0}, displayIters);
        model_m_fixed = vargplvmOptimiseModel(model_m_fixed, 0, 0, {0, iters}, displayIters);
        YpredBGPLVM_fixed_missing = modelPosteriorMeanVar(model_m_fixed, Xinterp);
    end
    if runGPLVM
        fprintf('\n\n# Optimising the GPLVM equivalent!\n\n');
        % fix all covars
        model_gplvm.vardist.covars = (model_gplvm.vardist.covars./model_gplvm.vardist.covars)*1e-10;
        model_gplvm.fixParamIndices = sort([trMeansInd trCovInd miCovInd]);
        vargplvmOptimiseModel(model_gplvm, 0, 0, {[ceil(initIters/4) ceil(initIters/4)], 0}, displayIters);
        model_gplvm = vargplvmOptimiseModel(model_gplvm, 0, 0, {0, iters}, displayIters);
        YpredGPLVM = modelPosteriorMeanVar(model_gplvm, Xinterp);
    end
end
%%
meanPred = repmat(mean(Yinterp), size(Yinterp,1), 1);

if ~isempty(resultsFile), outF = fopen(resultsFile, 'w'); else outF = ''; end
util_myFprintf(outF, '\n\n# ======== MSE ================\n');
try util_myFprintf(outF, '# Error mean                : %.3f\n', mean(abs(Yinterp(:) - meanPred(:)))); catch e, end
try util_myFprintf(outF, '# Error NN                  : %.3f\n', mean(abs(Yinterp(:) - YPredNN(:)))); catch e, end
try util_myFprintf(outF, '# Error NNx                 : %.3f\n', mean(abs(Yinterp(:) - YPredNNx(:)))); catch e, end
try util_myFprintf(outF, '# Error LinReg              : %.3f\n', mean(abs(Yinterp(:) - YPredLin(:)))); catch e, end
try util_myFprintf(outF, '# Error BGPLVM              : %.3f\n', mean(abs(Yinterp(:) - YpredBGPLVM(:)))); catch e, end
try util_myFprintf(outF, '# Error BGPLVM_missing      : %.3f\n', mean(abs(Yinterp(:) - YpredBGPLVM_missing(:)))); catch e, end
try util_myFprintf(outF, '# Error BGPLVM_fixed_missing: %.3f\n', mean(abs(Yinterp(:) - YpredBGPLVM_fixed_missing(:)))); catch e, end
try util_myFprintf(outF, '# Error GP                  : %.3f\n', mean(abs(Yinterp(:) - YpredGP(:)))); catch e, end
try util_myFprintf(outF, '# Error GPLVM               : %.3f\n', mean(abs(Yinterp(:) - YpredGPLVM(:)))); catch e, end
try util_myFprintf(outF, '# Error BGPLVM_init         : %.3f\n', mean(abs(Yinterp(:) - YpredBGPLVM_init(:)))); catch e, end


