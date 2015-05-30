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
% Modifications : Carl Henrik Ek, 2013
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
        fprintf('# Scaling output data to variance one...\n')
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
    if(length(options.scaleVal)==1)
        model.scale = repmat(options.scaleVal, 1, model.d);
    elseif(length(options.scaleVal)==model.d)
        model.scale = options.scaleVal;
    else
        error('Dimension mismatch in scaleVal');
    end
end

model.y = Y;
model.m = gpComputeM(model);

% Marg. likelihood is a sum ML = L + KL. Here we allow for
% ML = 2*((1-fw)*L + fw*KL), fw = options.KLweight. Default is 0.5, giving
% same weight to both terms.
model.KLweight = options.KLweight;
assert(model.KLweight >=0 && model.KLweight <=1);


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
model.learnBeta = 1;

% If the provided matrices are really big, the required quantities can be
% computed externally (e.g. block by block) and be provided here (we still
% need to add a switch to get that).
if model.d > limitDimensions && model.N < limitDimensions
    model.DgtN = 1; % D greater than N mode on.
    fprintf(1, '# The dataset has a large number of dimensions (%d)! Switching to "large D" mode!\n',model.d);
    
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
else
    % Trace(Y*Y) is a constant, so it can be calculated just once and stored
    % in memory to be used whenever needed.
    model.DgtN = 0;
    model.TrYY = sum(sum(model.m .* model.m));
end

model.date = date;

%%% Also, if there is a test dataset, Yts, then when we need to take
% model.m*my' (my=Yts_m(i,:)) in the pointLogLikeGradient, we can prepare in advance
% A = model.m*Yts_m' (NxN) and select A(:,i).
% Similarly, when I have m_new = [my;m] and then m_new*m_new', I can find that
% by taking m*m' (already computed as YY) and add one row on top: (m*my)'
% and one column on the left: m*my.
%
% %%%%%%%%%%



if isstruct(options.kern)
    model.kern = options.kern;
else
    model.kern = kernCreate(model.X, options.kern);
end



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
        model = vargplvmInitInducing(model, options);
        if isfield(options, 'beta')
            model.beta = options.beta;
        end
end


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

if isfield(options, 'learnInducing')
    model.learnInducing = options.learnInducing;
end

% Default field, controlling whether we are in initialising the var. distr.
% mode (by fixing SNR) or not.
model.initVardist = 0;

if isfield(options, 'fixVardist') && options.fixVardist
    model.fixParamIndices = 1:2*model.N*model.q;
    model.vardist.fixed = true;
end


%model = vargplvmCreate_extras(model, options, limitDimensions);

end

