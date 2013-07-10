function model = vargplvmInitDynamics(model,optionsDyn)
% VARGPLVMINITDYNAMICS Initialize the dynamics of a var-GPLVM model.
% FORMAT
% ARG optionsDyn : the VARGPLVM structure with the options used to
% create the vargplvm model..
%
%
% COPYRIGHT : Michalis K. Titsias, 2011
% COPYRIGHT : Neil D. Lawrence, 2011
% COPYRIGHT : Andreas C. Damianou, 2011
%
% SEEALSO : vargplvmOptionsDyn, vargplvmAddDynamics

% VARGPLVM



if ~isfield(model, 'dynamics') || isempty(model.dynamics)
    error(['model does not have field dynamics']);
end


%model.dynamics.kern.comp{1}.inverseWidth = 200./(((max(model.dynamics.t)-min(model.dynamics.t))).^2);
%params = vargplvmExtractParam(model);
%model = vargplvmExpandParam(model, params);


% Initialize barmu: first initialize X (mu) and then (see lines below)
% initialize barmu. The first 'if' statement checks if optionsDyn.initX is
% the initial X itself (e.g. calculated before in the demo) and there is no
% need to recalculate it.
if isfield(optionsDyn,'initX') && ~isstr(optionsDyn.initX)
    X = optionsDyn.initX;
else 
    if isfield(optionsDyn,'initX')
        initFunc = str2func([optionsDyn.initX 'Embed']);
    else
        initFunc=str2func('ppcaEmbed');
    end
    
    % If the model is in "D greater than N" mode, then we want initializations
    % to be performed with the original model.m
    if isfield(model, 'DgtN') && model.DgtN
        X = initFunc(model.mOrig, model.q);
    else
        X = initFunc(model.m, model.q);
    end
end


%---------------------- TEMP -------------------------------%
% The following are several tries to make the init. of the means better.
% It does that by "squeezing" the ends of the means to values closer to
% zero.

% X(1,:)=zeros(size(X(1,:)));
% X(end,:)=zeros(size(X(1,:)));

% X = 0.1*randn(size(X));
% X(1,:)=zeros(size(X(1,:)));
% X(end,:)=zeros(size(X(1,:)));



% m = floor((size(X,1))/2);
% p=max(max(abs(X))); p=p/10; %p = p^(1/m);
% l = logspace(p,0,m);
% l2 = logspace(p,0,m);
% l2 = l2(end:-1:1);
% l = [l l2]';
% mask = repmat(l,1,size(X,2));
% X = X.*(1./mask);
% X(1,:) = zeros(size(X(1,:)));
% X(end,:) = zeros(size(X(1,:)));
% X(end-1,:) = zeros(size(X(1,:)));

%save 'TEMPX.mat' 'X' %%TEMP

if optionsDyn.regularizeMeans
    model.dynamics.regularizeMeans = 1;
    m = floor((size(X,1))/7);  % 7
    p=ones(1,model.q);
    mm=max(X).^(1/2); %1/2 % the smallest this ratio(the exponent) is, the largest the effect
    for i=0:m-4
        X(end-i,:) = X(end-i,:) -  X(end-i,:).*p;
        p = p./mm;
    end
    p=ones(1,model.q);
    for i=1:m
        X(i,:) = X(i,:) -  X(i,:).*p;
        p = p./mm;
    end
end


%model.X = X;


% m = floor((size(X,1))/2);
% l = logspace(10,1,m);
% l2 = logspace(10,1,m);
% l2 = l2(end:-1:1);
% mask = repmat([l l2]',1,size(X,2));
% X2 = X-X.*mask;


%-------------------------------------------------------------------


vX = var(X);
noise = 0;%.01;
for q=1:model.q
    Lkt = chol(model.dynamics.Kt + noise*vX(q)*eye(model.N))';
    % barmu = inv(Kt + s*I)*X, so that  mu = Kt*barmu =  Kt*inv(Kt +
    % s*I)*X, which is jsut the GP prediction, is temporally smoothed version
    % of the PCA latent variables X (for the data Y)
    model.dynamics.vardist.means(:,q) = Lkt'\(Lkt\X(:,q));
end

if isfield(optionsDyn, 'vardistCovars') && ~isempty(optionsDyn.vardistCovars)
    if length(optionsDyn.vardistCovars) ==1
        model.dynamics.vardist.covars = 0.1*ones(model.N,model.q) + 0.001*randn(model.N,model.q);
        model.dynamics.vardist.covars(model.dynamics.vardist.covars<0.05) = 0.05;
        model.dynamics.vardist.covars = optionsDyn.vardistCovars * ones(size(model.dynamics.vardist.covars));
    elseif size(optionsDyn.vardistCovars) == size(model.dynamics.vardist.covars)
        model.dynamics.vardist.covars = optionsDyn.vardistCovars;
    end
else % Try to automatically find a satisfactory value for the initial reparametrized
     % covariances (lambda) so that the median of the real covariances is
     % as close as possible to an ideal value, eg 0.18
    bestVdc = 0.05; bestDiff = 1000;
    for vdc = [0.05:0.05:3]
        model.dynamics.vardist.covars = 0.1*ones(model.N,model.q) + 0.001*randn(model.N,model.q);
        model.dynamics.vardist.covars(model.dynamics.vardist.covars<0.05) = 0.05;
        model.dynamics.vardist.covars = vdc * ones(size(model.dynamics.vardist.covars));
        modelNew = vargplvmDynamicsUpdateStats(model);
        med = median(modelNew.vardist.covars(:));
        curDiff = abs(med - 0.18);
        if (curDiff) < bestDiff % "Ideal" median is 0.18
            bestVdc = vdc;
            bestDiff = curDiff;
        end
    end
    model.dynamics.vardist.covars = 0.1*ones(model.N,model.q) + 0.001*randn(model.N,model.q);
    model.dynamics.vardist.covars(model.dynamics.vardist.covars<0.05) = 0.05;
    model.dynamics.vardist.covars = bestVdc * ones(size(model.dynamics.vardist.covars));
    model = vargplvmDynamicsUpdateStats(model);
    fprintf('# Automatically calibrated lambda so that median(covars)=%.2f\n',median(model.vardist.covars(:)))
end

% smaller lengthscales for the base model
%model.kern.comp{1}.inputScales = 5./(((max(X)-min(X))).^2);
params = vargplvmExtractParam(model);
model = vargplvmExpandParam(model, params);


if isfield(optionsDyn,'X_u')
    model.X_u = optionsDyn.X_u;
else
    % inducing point need to initilize based on model.vardist.means
    if model.k <= model.N % model.N = size(mode.vardist.means,1) 
        if ~isfield(optionsDyn, 'labels')
            ind = randperm(model.N);
            ind = ind(1:model.k);
            model.X_u = model.X(ind, :);
        else
            % in the case that class labels are supplied, make sure that inducing inputs
            % from all classes are chosen
            [idcs, nSmpls] = class_samples( optionsDyn.labels, model.k );
            
            count = 1;
            midx = [];
            for inds = idcs
                ind   = inds{:};
                ind   = ind(randperm(numel(ind)));
                idx  = ind(1:nSmpls(count));
                
                % test that there is no overlap between index sets
                assert(isempty(intersect(midx, idx)));
                midx = [midx, idx];
                
                count = count+1;
            end
            model.X_u = model.X(midx,:);
        end
    else
        samplingInd=0; %% TEMP
        if samplingInd
            % !!! We could also try to sample all inducing points (more uniform
            % solution)
            % This only works if k<= 2*N
            model.X_u=zeros(model.k, model.q);
            ind = randperm(model.N);
            %ind = ind(1:model.N);
            model.X_u(1:model.N,:) = model.X(ind, :);
            
            % The remaining k-N points are sampled from the (k-N) first
            % distributions of the variational distribution (this could be done
            % randomly as well).
            dif=model.k-model.N;
            model.X_u(model.N+1:model.N+dif,:)=model.vardist.means(1:dif,:) + rand(size(model.vardist.means(1:dif,:))).*sqrt(model.vardist.covars(1:dif,:));  % Sampling from a Gaussian.
        else
            % !!! The following is not a good idea because identical initial
            % ind.points are difficult to optimise... better add some
            % additional noise.
            model.X_u=zeros(model.k, model.q);
            for i=1:model.k
                %ind=randi([1 size(model.vardist.means,1)]);
                % Some versions do not have randi... do it with rendperm
                % instead:
                ind=randperm(size(model.vardist.means,1));
                ind=ind(1);
                model.X_u(i,:) = model.vardist.means(ind,:)+0.001*randn(1,model.q);
            end
        end
    end
end

if isfield(optionsDyn, 'testReoptimise')
    model.dynamics.reoptimise = optionsDyn.testReoptimise;
end

model.dynamics.learnVariance = optionsDyn.learnVariance; % DEFAULT: 0

if isfield(optionsDyn, 'constrainType') && ~isempty(optionsDyn.constrainType)
    model.dynamics.constrainType = optionsDyn.constrainType;
end
    
params = vargplvmExtractParam(model);
model = vargplvmExpandParam(model, params);


med = median(model.vardist.covars(:));
minC = min(model.vardist.covars(:));
maxC = max(model.vardist.covars(:));
if med < 0.1 || med > 0.6
   warning('!!! [Min Median Max] value of variational covariances is [%1.2f %1.2f %1.2f].\n', minC, med, maxC);
else
    fprintf('# [Min Median Max] value of variational covariances is [%1.2f %1.2f %1.2f].\n', minC, med, maxC);
end
