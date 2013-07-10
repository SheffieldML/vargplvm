function g = svargplvmLogLikeGradients(model)

% SVARGPLVMLOGLIKEGRADIENTS Compute the gradients for the variational shared GPLVM.
% FORMAT
% DESC returns the gradients of the log likelihood with respect to the
% parameters of the GP-LVM model and with respect to the latent
% positions of the GP-LVM model.
% ARG model : the FGPLVM structure containing the parameters and
% the latent positions.
% RETURN g : the gradients of the latent positions (or the back
% constraint's parameters) and the parameters of the GP-LVM model.
%
% FORMAT
% DESC returns the gradients of the log likelihood with respect to the
% parameters of the GP-LVM model and with respect to the latent
% positions of the GP-LVM model in seperate matrices.
% ARG model : the FGPLVM structure containing the parameters and
% the latent positions.
% RETURN gX : the gradients of the latent positions (or the back
% constraint's parameters).
% RETURN gParam : gradients of the parameters of the GP-LVM model.
%
% COPYRIGHT : Andreas C. Damianou, 2011

% SEEALSO : svargplvmLogLikelihood, svargplvmCreate, modelLogLikeGradients

% VARGPLVM

try
    pool_open = matlabpool('size')>0;
catch e
    pool_open = 0;
end

if pool_open && (isfield(model,'parallel') && model.parallel)
    g=gPar(model);
else
    % Functions g1 and g2 should be equivalent for the static case, but
    % g1 might (?) be faster.
    if ~isfield(model, 'dynamics') || isempty(model.dynamics)
        g=g1(model);
        %g = g2(model);
    else
        if isfield(model.dynamics, 'seq') & ~isempty(model.dynamics.seq)
            g = gDyn(model); % Slower but memory efficient.
        else
            %g = gDyn(model); % Slower but memory efficient.
            g = gDynFast(model); % Faster, but memory consuming
        end
    end
end

%g = g1(model);
%g = g2(model);
end

%%%%%!!!!!! NOTE:
% g1 and g2 should be equivalent for the non-dynamics case. gPar is the
% equivalent of g1 for parallel computations. gDyn and gDynFast are
% suitable when there are dynamics, with the first focusing on memory
% efficiency and the second on speed.

function g = gPar(model)

% THE FOLLOWING WORKS ONLY WHEN THERE ARE NO DYNAMICS... %___ TODO TEMP
% Shared params
if isfield(model, 'dynamics') && isempty(model.dynamics) % TEMP
    error('The gradients for the dynamics case are not implemented correctly yet!'); % TEMP
else % ....
    gVarmeansKL = - model.vardist.means(:)';
    gVarcovsKL = 0.5 - 0.5*model.vardist.covars(:)';
end
model = svargplvmPropagateField(model, 'onlyLikelihood', 1);

%fprintf('# Derivs for KL is (should be zero): ');%%%TEMP

% Private params
% g = [[sum_m(gVar_only_likelihood)+gVar_onlyKL] g_1 g_2 ...]
% where sum_m is the sum over all models and g_m is the gradients for the
% non-shared parameters for model m
gShared = [gVarmeansKL gVarcovsKL];
modelTemp=model.comp;
parfor i=1:model.numModels
    gAll{i} = vargplvmLogLikeGradients(modelTemp{i});
end
g=[];
for i=1:model.numModels
    g_i = gAll{i};
    % Now add the derivatives for the shared parameters.
    if isfield(model, 'dynamics') & ~isempty(model.dynamics) % TEMP (doesn't work correctly for dynamics)
        gShared = gShared + g_i(1:model.dynamics.nParams); % TEMP (doesn't work correctly for dynamics)
        g_i = g_i((model.dynamics.nParams+1):end); % TEMP (doesn't work correctly for dynamics)
    else % else it's only the vardist. of the KL
        gShared = gShared + g_i(1:model.vardist.nParams);
        g_i = g_i((model.vardist.nParams+1):end);
    end
    g = [g g_i];
end
g = [gShared g];

end

function g = g1(model)


% THE FOLLOWING WORKS ONLY WHEN THERE ARE NO DYNAMICS... %___ TODO TEMP
% Shared params
if isfield(model, 'dynamics') && isempty(model.dynamics) % TEMP !!!! TODO (~isempty?)
    error('The gradients for the dynamics case are not implemented correctly yet!'); % TEMP
else % ....
    gVarmeansKL = - model.vardist.means(:)';
    gVarcovsKL = 0.5 - 0.5*model.vardist.covars(:)';
end
model = svargplvmPropagateField(model, 'onlyLikelihood', 1);

%fprintf('# Derivs for KL is (should be zero): ');%%%TEMP

% Private params
% g = [[sum_m(gVar_only_likelihood)+gVar_onlyKL] g_1 g_2 ...]
% where sum_m is the sum over all models and g_m is the gradients for the
% non-shared parameters for model m
g = [];
gShared = [gVarmeansKL gVarcovsKL];

for i=1:model.numModels
    g_i = vargplvmLogLikeGradients(model.comp{i});
    % Now add the derivatives for the shared parameters.
    if isfield(model, 'dynamics') & ~isempty(model.dynamics) % !! (doesn't work correctly for dynamics)
        gShared = gShared + g_i(1:model.dynamics.nParams); % !! (doesn't work correctly for dynamics)
        g_i = g_i((model.dynamics.nParams+1):end); % TEMP (doesn't work correctly for dynamics)
    else % else it's only the vardist. of the KL
        gShared = gShared + g_i(1:model.vardist.nParams);
        g_i = g_i((model.vardist.nParams+1):end);
    end
    g = [g g_i];
end
g = [gShared g];
end
% TEMP = vargplvmLogLikeGradients(model.comp{1}); %%%TEMP
% sum(sum(abs(gShared - TEMP(1:160)))) %%%TEMP
% g = [g 0*g_i]; %%%%%%%%%%% TEMP
%g = [-vargplvmLogLikeGradients(model.comp{1}) 0*g_i];


% NOT TESTED FOR THE DYNAMICS CASE
function g = g2(model)
g_1 = vargplvmLogLikeGradients(model.comp{1});
gShared = g_1(1:model.vardist.nParams);
g_1 = g_1(model.vardist.nParams+1:end);
model = svargplvmPropagateField(model, 'onlyLikelihood', 1, true);
g=[];
for i=2:model.numModels
    g_i = vargplvmLogLikeGradients(model.comp{i});
    % Now add the derivatives for the shared parameters.
    if isfield(model, 'dynamics') & ~isempty(model.dynamics)
        gShared = gShared + g_i(1:model.dynamics.nParams);
        g_i = g_i((model.dynamics.nParams+1):end);
    else % else it's only the vardist. of the KL
        gShared = gShared + g_i(1:model.vardist.nParams);
        g_i = g_i((model.vardist.nParams+1):end);
    end
    g = [g g_i];
end
g = [gShared g_1 g];
end



function g = gDynFast(modelAll)

gPrivAll = [];
gSharedCoeff = 0;
dynModel = modelAll.dynamics;
gSharedVar = zeros(1,dynModel.vardist.nParams);


for m=1:modelAll.numModels
    % This is similar to vargplvmLogLikeGradients but for every model
    
    model = modelAll.comp{m}; % current model
    
    
    % Include calculations for the KL term only once, for the first model.
    if m == 1
        includeKL = 1;
    else
        includeKL = 0;
    end
    
    
    
    % Likelihood terms (coefficients)
    [gK_uu, gPsi0, gPsi1, gPsi2, g_Lambda, gBeta, tmpV] = vargpCovGrads(model);
    
    if isfield(model, 'learnInducing')
        learnInducing = model.learnInducing;
    else
        learnInducing = true;
    end
    
    % Get (in three steps because the formula has three terms) the gradients of
    % the likelihood part w.r.t the data kernel parameters, variational means
    % and covariances (original ones). From the field model.vardist, only
    % vardist.means and vardist.covars and vardist.lantentDimension are used.
    [gKern1, gVarmeans1, gVarcovs1, gInd1] = kernVardistPsi1Gradient(model.kern, model.vardist, model.X_u, gPsi1', learnInducing);
    [gKern2, gVarmeans2, gVarcovs2, gInd2] = kernVardistPsi2Gradient(model.kern, model.vardist, model.X_u, gPsi2, learnInducing);
    [gKern0, gVarmeans0, gVarcovs0] = kernVardistPsi0Gradient(model.kern, model.vardist, gPsi0);
    gKern3 = kernGradient(model.kern, model.X_u, gK_uu);
    
    % At this point, gKern gVarmeansLik and gVarcovsLik have the derivatives for the
    % likelihood part. Sum all of them to obtain the final result.
    gKern = gKern0 + gKern1 + gKern2 + gKern3;
    gVarmeansLik = gVarmeans0 + gVarmeans1 + gVarmeans2;
    
    
    
    if strcmp(model.kern.type, 'rbfardjit')
        % different derivatives for the variance, which is super-numerically stable for
        % this particular kernel
        if model.learnSigmaf == 1
            gKern(1) = 0.5*model.d*( - model.k+ sum(sum(model.invLat.*model.invLat))/model.beta - model.beta*(model.Psi0-model.TrC)  )...
                + 0.5*tmpV;
            
            if ~isstruct(model.kern.transforms(1))
                fhandle = str2func([model.kern.transform(1) 'Transform']);
                gKern(1) = gKern(1).*fhandle(model.kern.variance, 'gradfact');
            else
                fhandle = str2func([model.kern.transforms(1).type 'Transform']);
                if ~isfield(model.kern.transforms(1), 'transformsettings')
                    gKern(1) = gKern(1).*fhandle(model.kern.variance, 'gradfact');
                else
                    gKern(1) = gKern(1).*fhandle(model.kern.variance, 'gradfact', model.kern.transforms(1).transformsettings);
                end
            end
        else
            gKern(1) = 0;
        end
    end
    
    %%% Compute Gradients with respect to X_u %%%
    gKX = kernGradX(model.kern, model.X_u, model.X_u);
    
    % The 2 accounts for the fact that covGrad is symmetric
    gKX = gKX*2;
    dgKX = kernDiagGradX(model.kern, model.X_u);
    for i = 1:model.k
        gKX(i, :, i) = dgKX(i, :);
    end
    
    if learnInducing
        % Allocate space for gX_u
        gX_u = zeros(model.k, model.q);
        % Compute portion associated with gK_u
        for i = 1:model.k
            for j = 1:model.q
                gX_u(i, j) = gKX(:, j, i)'*gK_uu(:, i);
            end
        end
        
        gInd = gInd1 + gInd2 + gX_u(:)';
    end
    
    % If the inducing points are fixed (tied to the latent points) then
    % X_u=K_t*dynamics.vardist.means and the derivatives w.r.t theta_t must be
    % amended with the appropriate partial derivatives. gInd must be passed,
    % in that case, as an argument to the function which calculates the
    % derivatives for the reparametrized quantities.
    if isfield(model, 'fixInducing') & model.fixInducing
        if learnInducing
            gIndRep = gInd;
        end
    else
        gIndRep=[];
    end
    
    gVarcovsLik = gVarcovs0 + gVarcovs1 + gVarcovs2;
    
    
    
    
    %----------------   This replaces vargpTimeDynamicsPriorReparamGrads
    
    % gVarmeansLik and gVarcovsLik are serialized into 1x(NxQ) vectors.
    % Convert them to NxQ matrices, i.e. fold them column-wise.
    gVarmeansLik = reshape(gVarmeansLik,dynModel.N,dynModel.q);
    gVarcovsLik = reshape(gVarcovsLik,dynModel.N,dynModel.q);
    
    if ~isempty(gIndRep)
        gIndRep = reshape(gIndRep, dynModel.N, dynModel.q); %%%%%
    end
    
    
    % The second term in the parenthesis, corresponds to the KL term. The
    % multiplier will set it to zero, if we need only the derivatives
    % corresponding only to the likelihood part.
    gVarmeans = dynModel.Kt * (gVarmeansLik - includeKL* dynModel.vardist.means);
    
    gVarcovs = zeros(dynModel.N, dynModel.q); % memory preallocation
    
    % The following quantities (in the loop) are needed to be calculated in
    % a single loop wrt q,
    if m==1
        for q=1:dynModel.q
            LambdaH_q = dynModel.vardist.covars(:,q).^0.5;
            Bt_q = eye(dynModel.N) + LambdaH_q*LambdaH_q'.*dynModel.Kt;
            
            % Invert Bt_q
            Lbt_q = jitChol(Bt_q)';
            G1 = Lbt_q \ diag(LambdaH_q);
            G = G1*dynModel.Kt;
            
            % Find Sq
            Sq{q} = dynModel.Kt - G'*G; %%%%%
            % Find the coefficient for the grad. wrt theta_t (params of Kt)
            G1T=G1';
            Bhat=G1T*G1;
            BhatKt=G1T*G;
            
            trGradKL{q} = -0.5*(BhatKt*Bhat + dynModel.vardist.means(:,q) * dynModel.vardist.means(:,q)');
            IBK = eye(dynModel.N) - BhatKt;
        end
    else
        for q=1:dynModel.q
            trGradKL{q} = 0;
        end
    end
    
    sumTrGradKL = 0;
    for q=1:dynModel.q
        % If includeKL is set to 0, then the KL part will be multiplied
        % with 0, thus being ignored.
        gVarcovs(:,q) = - (Sq{q} .* Sq{q}) * (gVarcovsLik(:,q) + includeKL * 0.5*dynModel.vardist.covars(:,q));
        
        diagVarcovs = repmat(gVarcovsLik(:,q)', dynModel.N,1);
        trGradKL{q} = trGradKL{q} + IBK .*  diagVarcovs * IBK';
        trGradKL{q} = trGradKL{q} + dynModel.vardist.means(:,q) * gVarmeansLik(:,q)';
        
        % In case gInd is empty then the inducing points are not reparametrized
        % (are not fixed to the variational means) we need not amend further the
        % derivative w.r.t theta_t, otherwise we have to do that.
        if ~isempty(gIndRep)
            trGradKL{q} = trGradKL{q} + dynModel.vardist.means(:,q) * gIndRep(:,q)';
        end
        
        sumTrGradKL = sumTrGradKL + trGradKL{q};
    end
    gSharedCoeff = gSharedCoeff + sumTrGradKL;
    
    
    
    
    % Serialize (unfold column-wise) gVarmeans and gVarcovs from NxQ matrices
    % to 1x(NxQ) vectors
    gVarcovs = gVarcovs(:)';
    gVarmeans = gVarmeans(:)';
    
    
    %-------------------------------------------------------------------
    
    
    
    % Variational variances are positive: Now that the final covariances
    % are obtained we amend with the partial derivatives due to the
    % exponential transformation to ensure positiveness.
    gVarcovs = (gVarcovs(:).*model.dynamics.vardist.covars(:))';
    
    if isfield(model, 'fixInducing') && model.fixInducing && learnInducing
        % If there are dynamics the derivative must further be amended with the
        % partial deriv. due to the mean reparametrization.
        if isfield(model, 'dynamics') && ~isempty(model.dynamics)
            gInd = reshape(gInd,model.k,model.q);
            %gInd = gInd' * model.dynamics.Kt;
            gInd =  model.dynamics.Kt * gInd;
            gInd = gInd(:)';
        end
        %gVarmeans(model.inducingIndices, :) = gVarmeans(model.inducingIndices,
        %:) + gInd; % This should work AFTER reshaping the matrices...but here
        %we use all the indices anyway.
        gVarmeans = gVarmeans + gInd;
        gInd = []; % Inducing points are not free variables anymore, they dont have derivatives on their own.
    end
    
    gVar = [gVarmeans gVarcovs];
    
    if isfield(model.vardist,'paramGroups')
        gVar = gVar*model.vardist.paramGroups;
    end
    
    % In case we are in the phase where the vardistr. is initialised (see above
    % for the variance of the kernel), beta is kept fixed. For backwards
    % compatibility this can be controlled either with the learnBeta field or
    % with the initVardist field. The later overrides the first.
    if isfield(model, 'learnBeta') && model.learnBeta
        gBetaFinal = gBeta;
    else
        gBetaFinal = 0*gBeta;
    end
    if isfield(model, 'initVardist')
        if model.initVardist == 1
            gBetaFinal = 0*gBeta;
        else
            gBetaFinal = gBeta;
        end
    end
    
    gSharedVar = gSharedVar + gVar;
    
    if ~learnInducing
        gInd = [];
    end
    
    
    % At this point, gDynKern will be [] if there are no dynamics.
    gPrivAll = [gPrivAll gInd gKern gBetaFinal];
end

gDynKern = kernGradient(dynModel.kern, dynModel.t, gSharedCoeff);


%%%% ...... was commented...
if isfield(model, 'dynamics') && ~isempty(model.dynamics)
    if strcmp(model.dynamics.kern.comp{1}.type,'rbf') || strcmp(model.dynamics.kern.comp{1}.type,'matern32') || strcmp(model.dynamics.kern.comp{1}.type,'rbfperiodic') || strcmp(model.dynamics.kern.comp{1}.type,'rbfperiodic2')
        if ~isfield(model.dynamics, 'learnVariance') || ~model.dynamics.learnVariance
            gDynKern(2) = 0;
        end
    end
    
    %___NEW: assume that the second rbf/matern etc kernel is last in the
    %compound kernel
    %if numel(model.dynamics.kern.comp) > 3
    if isfield(model.dynamics, 'learnSecondVariance') && ~model.dynamics.learnSecondVariance   %%%%% NEW
        gDynKern(end) = 0;
    end
    %end
    %___
end


%---- NEW 2012: This is to fix selectively some of the kernel's parameters
if ~isfield(model.dynamics, 'learnVariance') || ~model.dynamics.learnVariance
    % The field model.dynamics.fixedVariance must have the
    % indexes of the gradient vector that correspond to
    % variances of kernels of the matern class and we wish to
    % not learn them.
    if isfield(model.dynamics, 'fixedKernVariance') && ~isempty(model.dynamics.fixedKernVariance)
        gDynKern(model.dynamics.fixedKernVariance) = 0;
    end
end
%----


g = [gSharedVar gDynKern gPrivAll];

end

% The bound is written as:
% ll = KL + ll1 + ll2 + ...
% For the derivatives, it's a bit trickier, since the functions we have
% from vargplvm calculate the derivatives for ll+KL in a mixed way (because
% the derivatives w.r.t mu_bar and lambda exist ALSO in the likelihood
% parts, due to reparametrising mu and S with Kt. So, we split as follows:
% [gVarAll gDynKern gPriv1 gPriv2 ... ]
% TODO: for varcovs, perform (Sq.*Sq)*[large sum]
%   and for varmeans, perform Kt * [large sum
% TODO: add switch so that one model can be dynamical, another can be
%   static
function g = gDyn(modelAll)

gPrivAll = [];
gSharedCoeff = 0;
dynModel = modelAll.dynamics;
gSharedVar = zeros(1,dynModel.vardist.nParams);

if isfield(dynModel, 'seq') & ~isempty(dynModel.seq)
    seqStart=1;
    seq = dynModel.seq;
    for i=1:length(dynModel.seq)
        seqEnd = seq(i);
        sumTrGradKL{i} = zeros(seqEnd-seqStart+1, seqEnd-seqStart+1);
        seqStart = seqEnd+1;
    end
end

for m=1:modelAll.numModels
    model = modelAll.comp{m}; % current model
    
    % Include calculations for the KL term only once, for the first model.
    if m == 1
        includeKL = 1;
    else
        includeKL = 0;
    end
    
    % Likelihood terms (coefficients)
    [gK_uu, gPsi0, gPsi1, gPsi2, g_Lambda, gBeta, tmpV] = vargpCovGrads(model);
    
    % Get (in three steps because the formula has three terms) the gradients of
    % the likelihood part w.r.t the data kernel parameters, variational means
    % and covariances (original ones). From the field model.vardist, only
    % vardist.means and vardist.covars and vardist.lantentDimension are used.
    [gKern1, gVarmeans1, gVarcovs1, gInd1] = kernVardistPsi1Gradient(model.kern, model.vardist, model.X_u, gPsi1');
    [gKern2, gVarmeans2, gVarcovs2, gInd2] = kernVardistPsi2Gradient(model.kern, model.vardist, model.X_u, gPsi2);
    [gKern0, gVarmeans0, gVarcovs0] = kernVardistPsi0Gradient(model.kern, model.vardist, gPsi0);
    gKern3 = kernGradient(model.kern, model.X_u, gK_uu);
    
    % At this point, gKern gVarmeansLik and gVarcovsLik have the derivatives for the
    % likelihood part. Sum all of them to obtain the final result.
    gKern = gKern0 + gKern1 + gKern2 + gKern3;
    gVarmeansLik = gVarmeans0 + gVarmeans1 + gVarmeans2;
    
    if strcmp(model.kern.type, 'rbfardjit')
        % different derivatives for the variance, which is super-numerically stable for
        % this particular kernel
        if model.learnSigmaf == 1
            gKern(1) = 0.5*model.d*( - model.k+ sum(sum(model.invLat.*model.invLat))/model.beta - model.beta*(model.Psi0-model.TrC)  )...
                + 0.5*tmpV;
            
            if ~isstruct(model.kern.transforms(1))
                fhandle = str2func([model.kern.transform(1) 'Transform']);
                gKern(1) = gKern(1).*fhandle(model.kern.variance, 'gradfact');
            else
                fhandle = str2func([model.kern.transforms(1).type 'Transform']);
                if ~isfield(model.kern.transforms(1), 'transformsettings')
                    gKern(1) = gKern(1).*fhandle(model.kern.variance, 'gradfact');
                else
                    gKern(1) = gKern(1).*fhandle(model.kern.variance, 'gradfact', model.kern.transforms(1).transformsettings);
                end
            end
        else
            gKern(1) = 0;
        end
    end
    
    %%% Compute Gradients with respect to X_u %%%
    gKX = kernGradX(model.kern, model.X_u, model.X_u);
    
    % The 2 accounts for the fact that covGrad is symmetric
    gKX = gKX*2;
    dgKX = kernDiagGradX(model.kern, model.X_u);
    for i = 1:model.k
        gKX(i, :, i) = dgKX(i, :);
    end
    
    % Allocate space for gX_u
    gX_u = zeros(model.k, model.q);
    % Compute portion associated with gK_u
    for i = 1:model.k
        for j = 1:model.q
            gX_u(i, j) = gKX(:, j, i)'*gK_uu(:, i);
        end
    end
    
    gInd = gInd1 + gInd2 + gX_u(:)';
    
    % If the inducing points are fixed (tied to the latent points) then
    % X_u=K_t*dynamics.vardist.means and the derivatives w.r.t theta_t must be
    % amended with the appropriate partial derivatives. gInd must be passed,
    % in that case, as an argument to the function which calculates the
    % derivatives for the reparametrized quantities.
    if isfield(model, 'fixInducing') & model.fixInducing
        gIndRep = gInd;
    else
        gIndRep=[];
    end
    gVarcovsLik = gVarcovs0 + gVarcovs1 + gVarcovs2;
    
    
    
    %----------------   This replaces vargpTimeDynamicsPriorReparamGrads
    
    % gVarmeansLik and gVarcovsLik are serialized into 1x(NxQ) vectors.
    % Convert them to NxQ matrices, i.e. fold them column-wise.
    gVarmeansLik = reshape(gVarmeansLik,dynModel.N,dynModel.q);
    gVarcovsLik = reshape(gVarcovsLik,dynModel.N,dynModel.q);
    
    if ~isempty(gIndRep)
        gIndRep = reshape(gIndRep, dynModel.N, dynModel.q); %%%%%
    end
    
    % The second term in the parenthesis, corresponds to the KL term. The
    % multiplier will set it to zero, if we need only the derivatives
    % corresponding only to the likelihood part.
    gVarmeans = dynModel.Kt * (gVarmeansLik - includeKL* dynModel.vardist.means);
    gVarcovs = zeros(dynModel.N, dynModel.q); % memory preallocation
    
    if isfield(dynModel, 'seq') & ~isempty(dynModel.seq)
        for q=1:dynModel.q
            LambdaH_q = dynModel.vardist.covars(:,q).^0.5;
            
            % The calculations are performed for each sequence independently
            % because correlations between sequences are not captured and, thus,
            % most of the matrices involved are block diagonal, so that in each of
            % the following loops only one block is considered.
            seqStart=1;
            for i=1:length(dynModel.seq)
                seqEnd = seq(i);
                Bt_q = eye(seqEnd-seqStart+1) + LambdaH_q(seqStart:seqEnd,1)*LambdaH_q(seqStart:seqEnd,1)'.*dynModel.Kt(seqStart:seqEnd,seqStart:seqEnd);
                
                % Invert Bt_q
                Lbt_q = jitChol(Bt_q)';
                G1 = Lbt_q \ diag(LambdaH_q(seqStart:seqEnd,1));
                G = G1*dynModel.Kt(seqStart:seqEnd, seqStart:seqEnd);
                
                % Find Sq
                Sq = dynModel.Kt(seqStart:seqEnd, seqStart:seqEnd) - G'*G;
                
                gVarcovs(seqStart:seqEnd,q) = - (Sq .* Sq) * (gVarcovsLik(seqStart:seqEnd,q) + ...
                    includeKL* 0.5*dynModel.vardist.covars(seqStart:seqEnd,q));
                
                % Find the coefficient for the grad. wrt theta_t (params of Kt)
                G1T=G1';
                Bhat=G1T*G1;
                BhatKt=G1T*G;
                
                if includeKL
                    trGradKL = -0.5*(BhatKt*Bhat + dynModel.vardist.means(seqStart:seqEnd,q) * dynModel.vardist.means(seqStart:seqEnd,q)');
                else
                    trGradKL = 0;
                end
                IBK = eye(seqEnd-seqStart+1) - BhatKt;
                diagVarcovs = repmat(gVarcovsLik(seqStart:seqEnd,q)', seqEnd-seqStart+1,1);
                trGradKL = trGradKL + IBK .*  diagVarcovs * IBK';
                trGradKL = trGradKL + dynModel.vardist.means(seqStart:seqEnd,q) * gVarmeansLik(seqStart:seqEnd,q)';
                
                % In case gInd is empty then the inducing points are not reparametrized
                % (are not fixed to the variational means) we need not amend further the
                % derivative w.r.t theta_t.
                if ~isempty(gIndRep)
                    trGradKL = trGradKL + dynModel.vardist.means(seqStart:seqEnd,q) * gIndRep(seqStart:seqEnd,q)';
                end
                sumTrGradKL{i} = sumTrGradKL{i} + trGradKL;
                seqStart = seqEnd+1;
            end
            %sumTrGradKL = sumTrGradKL + trGradKL;
        end
    else
        sumTrGradKL = 0;
        for q=1:dynModel.q
            LambdaH_q = dynModel.vardist.covars(:,q).^0.5;
            Bt_q = eye(dynModel.N) + LambdaH_q*LambdaH_q'.*dynModel.Kt;
            
            % Invert Bt_q
            Lbt_q = jitChol(Bt_q)';
            G1 = Lbt_q \ diag(LambdaH_q);
            G = G1*dynModel.Kt;
            
            % Find Sq
            Sq = dynModel.Kt - G'*G;
            
            % If onlyLikelihood is set to 1, then the KL part will be multiplied
            % with zero, thus being ignored.
            gVarcovs(:,q) = - (Sq .* Sq) * (gVarcovsLik(:,q) + includeKL * 0.5*dynModel.vardist.covars(:,q));
            
            % Find the coefficient for the grad. wrt theta_t (params of Kt)
            G1T=G1';
            Bhat=G1T*G1;
            BhatKt=G1T*G;
            
            if includeKL
                trGradKL = -0.5*(BhatKt*Bhat + dynModel.vardist.means(:,q) * dynModel.vardist.means(:,q)');
            else
                trGradKL = 0;
            end
            
            IBK = eye(dynModel.N) - BhatKt;
            diagVarcovs = repmat(gVarcovsLik(:,q)', dynModel.N,1);
            trGradKL = trGradKL + IBK .*  diagVarcovs * IBK';
            trGradKL = trGradKL + dynModel.vardist.means(:,q) * gVarmeansLik(:,q)';
            
            % In case gInd is empty then the inducing points are not reparametrized
            % (are not fixed to the variational means) we need not amend further the
            % derivative w.r.t theta_t, otherwise we have to do that.
            if ~isempty(gIndRep)
                trGradKL = trGradKL + dynModel.vardist.means(:,q) * gIndRep(:,q)';
            end
            sumTrGradKL = sumTrGradKL + trGradKL;
        end
        gSharedCoeff = gSharedCoeff + sumTrGradKL;
    end
    % Serialize (unfold column-wise) gVarmeans and gVarcovs from NxQ matrices
    % to 1x(NxQ) vectors
    gVarcovs = gVarcovs(:)';
    gVarmeans = gVarmeans(:)';
    %-------------------------------------------------------------------
    
    
    % Variational variances are positive: Now that the final covariances
    % are obtained we amend with the partial derivatives due to the
    % exponential transformation to ensure positiveness.
    gVarcovs = (gVarcovs(:).*model.dynamics.vardist.covars(:))';
    
    if isfield(model, 'fixInducing') & model.fixInducing
        % If there are dynamics the derivative must further be amended with the
        % partial deriv. due to the mean reparametrization.
        if isfield(model, 'dynamics') && ~isempty(model.dynamics)
            gInd = reshape(gInd,model.k,model.q);
            %gInd = gInd' * model.dynamics.Kt;
            gInd =  model.dynamics.Kt * gInd;
            gInd = gInd(:)';
        end
        %gVarmeans(model.inducingIndices, :) = gVarmeans(model.inducingIndices,
        %:) + gInd; % This should work AFTER reshaping the matrices...but here
        %we use all the indices anyway.
        gVarmeans = gVarmeans + gInd;
        gInd = []; % Inducing points are not free variables anymore, they dont have derivatives on their own.
    end
    gVar = [gVarmeans gVarcovs];
    
    if isfield(model.vardist,'paramGroups')
        gVar = gVar*model.vardist.paramGroups;
    end
    
    if isfield(model, 'dynamics') && ~isempty(model.dynamics)
        if strcmp(model.dynamics.kern.comp{1}.type,'rbf') || strcmp(model.dynamics.kern.comp{1}.type,'matern32') || strcmp(model.dynamics.kern.comp{1}.type,'rbfperiodic') || strcmp(model.dynamics.kern.comp{1}.type,'rbfperiodic2')
            if ~isfield(model.dynamics, 'learnVariance') || ~model.dynamics.learnVariance
                gDynKern(2) = 0;
            end
        end
        %___NEW: assume that the second rbf/matern etc kernel is last in the
        %compound kernel
        %if numel(model.dynamics.kern.comp) > 3
        if isfield(model.dynamics, 'learnSecondVariance') && ~model.dynamics.learnSecondVariance   %%%%% NEW
            gDynKern(end) = 0;
        end
        %end
        %___
    end
    
    % In case we are in the phase where the vardistr. is initialised (see above
    % for the variance of the kernel), beta is kept fixed. For backwards
    % compatibility this can be controlled either with the learnBeta field or
    % with the initVardist field. The later overrides the first.
    if isfield(model, 'learnBeta') && model.learnBeta
        gBetaFinal = gBeta;
    else
        gBetaFinal = 0*gBeta;
    end
    if isfield(model, 'initVardist')
        if model.initVardist == 1
            gBetaFinal = 0*gBeta;
        else
            gBetaFinal = gBeta;
        end
    end
    
    gSharedVar = gSharedVar + gVar;
    % At this point, gDynKern will be [] if there are no dynamics.
    gPrivAll = [gPrivAll gInd gKern gBetaFinal];
end

if isfield(dynModel, 'seq') & ~isempty(dynModel.seq)
    seqStart=1;
    gDynKern = size(dynModel.kern.nParams,1);
    for i=1:length(dynModel.seq)
        seqEnd = seq(i);
        gDynKern = gDynKern + kernGradient(dynModel.kern, dynModel.t(seqStart:seqEnd), sumTrGradKL{i});
        seqStart = seqEnd+1;
    end
else
    gDynKern = kernGradient(dynModel.kern, dynModel.t, gSharedCoeff);
end
g = [gSharedVar gDynKern gPrivAll];

end



function [gK_uu, gPsi0, gPsi1, gPsi2, g_Lambda, gBeta, tmpV] = vargpCovGrads(model)

gPsi1 = model.beta * model.m * model.B';
gPsi1 = gPsi1'; % because it is passed to "kernVardistPsi1Gradient" as gPsi1'...

gPsi2 = (model.beta/2) * model.T1;

gPsi0 = -0.5 * model.beta * model.d;

gK_uu = 0.5 * (model.T1 - (model.beta * model.d) * model.invLmT * model.C * model.invLm);

sigm = 1/model.beta; % beta^-1

PLm = model.invLatT*model.P;
tmpV = sum(sum(PLm.*PLm));
gBeta = 0.5*(model.d*(model.TrC + (model.N-model.k)*sigm -model.Psi0) ...
    - model.TrYY + model.TrPP ...
    + (1/(model.beta^2)) * model.d * sum(sum(model.invLat.*model.invLat)) + sigm*tmpV);


if ~isstruct(model.betaTransform)
    fhandle = str2func([model.betaTransform 'Transform']);
    gBeta = gBeta*fhandle(model.beta, 'gradfact');
else
    fhandle = str2func([model.betaTransform.type 'Transform']);
    gBeta = gBeta*fhandle(model.beta, 'gradfact', model.betaTransform.transformsettings);
end

g_Lambda = repmat(-0.5*model.beta*model.d, 1, model.N);
end
