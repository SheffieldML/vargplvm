function model = vargplvmCreate(q, d, Y, options)

% VARGPLVMCREATE Create a GPLVM model with inducing variables.
% FORMAT
% DESC creates a GP-LVM model with the possibility of using
% inducing variables to speed up computation.
% ARG q : dimensionality of latent space.
% ARG d : dimensionality of data space.
% ARG Y : the data to be modelled in design matrix format (as many
% rows as there are data points).
% ARG options : options structure as returned from
% VARGPLVMOPTIONS. This structure determines the type of
% approximations to be used (if any).
% ARG enableDgtN: if set to true, the model will be created in D >> N
% mode (if this is indeed the case).
% RETURN model : the GP-LVM model.
%
% COPYRIGHT : Michalis K. Titsias and Neil D. Lawrence, 2009-2011
% Modifications : Andreas C. Damianou, 2010-2011
%
% SEEALSO : vargplvmOptions

% VARGPLVM

if size(Y, 2) ~= d
    error(['Input matrix Y does not have dimension ' num2str(d)]);
end


if isfield(options,'enableDgtN')
	enableDgtN = options.enableDgtN;
else
	enableDgtN = true; % default behaviour is to enable D >> N mode
end

% Datasets with dimensions larger than this number will run the
% code which is more efficient for large D. This will happen only if N is
% smaller than D, though.
if enableDgtN
    limitDimensions = 2500; % Default: 5000
else
    limitDimensions = 1e+10;
end

model.type = 'vargplvm';
model.approx = options.approx;

model.learnScales = options.learnScales;
%model.scaleTransform = optimiDefaultConstraint('positive');

if isfield(options, 'optimiseBeta')
    model.optimiseBeta = options.optimiseBeta;
end

if isfield(options, 'betaTransform') && isstruct( options.betaTransform )
    model.betaTransform = options.betaTransform;
else
    model.betaTransform =  optimiDefaultConstraint('positive');
end

model.q = q;
model.d = size(Y, 2);
model.N = size(Y, 1);

model.optimiser = options.optimiser;
model.bias = mean(Y);
model.scale = ones(1, model.d);

if(isfield(options,'scale2var1'))
    if(options.scale2var1)
        model.scale = std(Y);
        model.scale(find(model.scale==0)) = 1;
        if(model.learnScales)
            warning('Both learn scales and scale2var1 set for GP');
        end
        if(isfield(options, 'scaleVal'))
            warning('Both scale2var1 and scaleVal set for GP');
        end
    end
end
if(isfield(options, 'scaleVal'))
    model.scale = repmat(options.scaleVal, 1, model.d);
end

model.y = Y;
model.m = gpComputeM(model);


%if options.computeS
%  model.S = model.m*model.m';
%  if ~strcmp(model.approx, 'ftc')
%    error('If compute S is set, approximation type must be ''ftc''')
%  end
%end



if isstr(options.initX)
    %%% The following should eventually be uncommented so as to initialize
    %%% in the dual space. This is much more efficient for large D.
    %if model.d > limitDimensions && model.N < limitDimensions
    %      [X vals] = svd(model.m*model.m');
    %      X = X(:,1:model.q); %%%%%
    %else
    initFunc = str2func([options.initX 'Embed']);
    X = initFunc(model.m, q);
    %end
else
    if size(options.initX, 1) == size(Y, 1) ...
            & size(options.initX, 2) == q
        X = options.initX;
    else
        error('options.initX not in recognisable form.');
    end
end

model.X = X;


model.learnBeta = 1; % Newly added: deafault value for learnBeta.

% If the provided matrices are really big, the required quantities can be
% computed externally (e.g. block by block) and be provided here (we still
% need to add a switch to get that).
if model.d > limitDimensions && model.N < limitDimensions
    model.DgtN = 1; % D greater than N mode on.
    fprintf(1, '# The dataset has a large number of dimensions (%d)! Switching to "large D" mode!\n',model.d);

    % If we have test data, we can prepare some multiplications in
    % advance, obtain NxN matrices and then never store Y.

    %{
    % Also, check that the following trick must only be done, in case there
    % is Yts, when length(indexPresent) > limitDimensions
    if(isfield(options,'Ytest'))
        % The same indices must be missing from all test points.
        indexMissing = find(isnan(options.Ytest(1,:)));
        indexPresent = setdiff(1:model.d, indexMissing);
        % Scale test data exactly as the training ones.
        modelTmp = model;
        modelTmp.y = options.Ytest; %%% Temporary assignment
        modelTmp.y = model.y(:,indexPresent);
        modelTmp.m = model.m(:,indexPresent);
        modelTmp.d = length(indexPresent);
        my = gpComputeM(model);
        modelTmp.y = []; % Never used
        my = my(:,indexPresent);


        mExt = [modelTmp.m; my];
        mExtmExtT = mExt * mExt';
        model.mExt = chol(mExtmExtT, 'lower'); % NxN
        model.TrMMExt = sum(diag(mExtmExtT)); % scalar
        model.mmyT = model.mExt * my'; % NxN
    end
    %}

    % yyT_denorm=model.y * model.y';
    % Keep the original m. It is needed for predictions.
    model.mOrig = model.m;

    % The following will eventually be uncommented.
    %model.y = []; % Never used.
    YYT = model.m * model.m'; % NxN
    % Replace data with the cholesky of Y*Y'.Same effect, since Y only appears as Y*Y'.
    %%% model.m = chol(YYT, 'lower');  %%% Put a switch here!!!!
    [U S V]=svd(YYT);
    model.m=U*sqrt(abs(S));

    model.TrYY = sum(diag(YYT)); % scalar

    %{
%     % Scale and bias
%     d2=size(model.m,2);
%     model.scale = ones(1, d2);
%
%     [U S V]=svd(yyT_denorm);
%     y2=U*sqrt(abs(S));
%     model.bias = mean(y2);
%     if(isfield(options,'scale2var1'))
%         if(options.scale2var1)
%             model.scale = std(y2);
%             model.scale(find(model.scale==0)) = 1;
%             if(model.learnScales)
%                 warning('Both learn scales and scale2var1 set for GP');
%             end
%             if(isfield(options, 'scaleVal'))numel(cidx{count});
%                 warning('Both scale2var1 and scaleVal set for GP');
%             end
%         end
%     end
%     if(isfield(options, 'scaleVal'))
%         model.scale = repmat(options.scaleVal, 1, d2);
%     end
    %}
else
    % Trace(Y*Y) is a constant, so it can be calculated just once and stored
    % in memory to be used whenever needed.
    model.DgtN = 0;
    model.TrYY = sum(sum(model.m .* model.m));
end

model.date = date;%%%%

%%% Also, if there is a test dataset, Yts, then when we need to take
% model.m*my' (my=Yts_m(i,:)) in the pointLogLikeGradient, we can prepare in advance
% A = model.m*Yts_m' (NxN) and select A(:,i).
% Similarly, when I have m_new = [my;m] and then m_new*m_new', I can find that
% by taking m*m' (already computed as YY) and add one row on top: (m*my)'
% and one column on the left: m*my.
%
% %%%%%%%%%%

%%%% _NEW



%model.isMissingData = options.isMissingData;
%if model.isMissingData
%  for i = 1:model.d
%    model.indexPresent{i} = find(~isnan(y(:, i)));
% end
%end
%model.isSpherical = options.isSpherical;


if isstruct(options.kern)
    model.kern = options.kern;
else
    model.kern = kernCreate(model.X, options.kern);
end

%{
%if isfield(options, 'noise')
%  if isstruct(options.noise)
%    model.noise = options.noise;
%  else
%    model.noise = noiseCreate(options.noise, y);
%  end
%
%  % Set up noise model gradient storage.
%  model.nu = zeros(size(y));
%  model.g = zeros(size(y));
%  model.gamma = zeros(size(y));
%
%  % Initate noise model
%  model.noise = noiseCreate(noiseType, y);numel(cidx{count});
%
%  % Set up storage for the expectations
%  model.expectations.f = model.y;
%  model.expectations.ff = ones(size(model.y));
%  model.expectations.fBar =ones(size(model.y));
%  model.expectations.fBarfBar = ones(numData, ...
%                                     numData, ...
%                                     size(model.y, 2));
%end
%}

% check if parameters are to be optimised in model space (and constraints are handled
% by optimiser)
if isfield(options, 'notransform') && options.notransform == true
    % store notransform option in model:
    % this is needed when computing gradients
    model.notransform = true;
    model.vardist = vardistCreate(X, q, 'gaussian', 'identity');
else
    model.vardist = vardistCreate(X, q, 'gaussian'); %%%
end

switch options.approx
    case {'dtcvar'}
        % Sub-sample inducing variables.
        model.k = options.numActive;
        model.fixInducing = options.fixInducing;
        if options.fixInducing
            if length(options.fixIndices)~=options.numActive
                %error(['Length of indices for fixed inducing variables must ' ...
                %    'match number of inducing variables']);
                warning(['Length of indices for fixed inducing variables was NOT' ...
                    ' equal to the number of inducing variables. This was automatically' ...
                    ' fixed by custom fixIndices.']);
                options.fixIndices = 1:options.numActive;
            end
            model.X_u = model.X(options.fixIndices, :);
            model.inducingIndices = options.fixIndices;
        else
            %%%NEW_: make it work even if k>N
            if model.k <= model.N
                if ~isfield(options, 'labels')
                    if ~isfield(options, 'initX_u')
                        ind = randperm(model.N);
                        ind = ind(1:model.k);
                        model.X_u = model.X(ind, :);
                    else
                        model.X_u = options.initX_u;
                    end
                else
                    % in the case that class labels are supplied, make sure that inducing inputs
                    % from all classes are chosen
                    [idcs, nSmpls] = class_samples( options.labels, model.k );
                    
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
                warning('k > N !')
                % TODO: sample from the variational distr. (this should probably go
                % to the dynamics as well because the vardist. changes in the initialization for the dynamics.

                %!!! The following code needs some more testing!
                samplingInd=0; %% TEMP
                if samplingInd
                    % This only works if k<= 2*N
                    model.X_u=zeros(model.k, model.q);
                    ind = randperm(model.N);
                    %ind = ind(1:model.N);
                    model.X_u(1:model.N,:) = model.X(ind, :);

                    % The remaining k-N points are sampled from the (k-N) first
                    % distributions of the variational distribution (this could be done
                    % randomly as well).
                    dif=model.k-model.N;
                    model.X_u(model.N+1:model.N+dif,:)= model.vardist.means(1:dif,:) + rand(size(model.vardist.means(1:dif,:))).*sqrt(model.vardist.covars(1:dif,:));  % Sampling from a Gaussian.
                else
                    model.X_u=zeros(model.k, model.q);
                    for i=1:model.k
                        %ind=randi([1 size(model.vardist.means,1)]);
                        % Some versions do not have randi... do it with rendperm
                        % instead: 
                        % ceil(size(model.vardist.means,1).*rand) % alternative
                        ind=randperm(size(model.vardist.means,1));
                        ind=ind(1);
                        model.X_u(i,:) = model.vardist.means(ind,:);
                    end
                end
            end
            %%%_NEW
        end
        if isfield(options, 'beta')
            model.beta = options.beta;
        end
end
%{
%if model.k>model.N
%  error('Number of active points cannot be greater than number of data.')
%end
%if strcmp(model.approx, 'pitc')
%  numBlocks = ceil(model.N/model.k);
%  numPerBlock = ceil(model.N/numBlocks);
%  startVal = 1;
%  endVal = model.k;
%  model.blockEnd = zeros(1, numBlocks);
%  for i = 1:numBlocks
%    model.blockEnd(i) = endVal;
%    endVal = numPerBlock + endVal;
%    if endVal>model.N
%      endVal = model.N;
%    end
%  end
%end
%}

if isstruct(options.prior)
    model.prior = options.prior;
else
    if ~isempty(options.prior)
        model.prior = priorCreate(options.prior);
    end
end

if isfield(options, 'notransform') && options.notransform == true
    model.prior.transforms.type = 'identity';    
end

%model.vardist = vardistCreate(X, q, 'gaussian');

if isfield(options, 'tieParam') & ~isempty(options.tieParam)
    %
    if strcmp(options.tieParam,'free')
        % paramsList =
    else
        startVal = model.vardist.latentDimension*model.vardist.numData + 1;
        endVal = model.vardist.latentDimension*model.vardist.numData;
        for q=1:model.vardist.latentDimension
            endVal = endVal + model.vardist.numData;
            index = startVal:endVal;
            paramsList{q} = index;
            startVal = endVal + 1;
        end
        model.vardist = modelTieParam(model.vardist, paramsList);
    end
    %
end

%{
%  if isstruct(options.back)

%if isstruct(options.inducingPrior)
%  model.inducingPrior = options.inducingPrior;
%else
%  if ~isempty(options.inducingPrior)
%    model.inducingPrior = priorCreate(options.inducingPrior);
%  end
%end

%if isfield(options, 'back') & ~isempty(options.back)
%  if isstruct(options.back)
%    model.back = options.back;
%  else
%    if ~isempty(options.back)
%      model.back = modelCreate(options.back, model.d, model.q, options.backOptions);
%    end
%  end
%  if options.optimiseInitBack
%    % Match back model to initialisation.
%    model.back = mappingOptimise(model.back, model.y, model.X);
%  end
%  % Now update latent positions with the back constraints output.
%  model.X = modelOut(model.back, model.y);
%else
%  model.back = [];
%end

model.constraints = {};

model.dynamics = [];

initParams = vargplvmExtractParam(model);
model.numParams = length(initParams);

% This forces kernel computation.
model = vargplvmExpandParam(model, initParams);
%}

end

