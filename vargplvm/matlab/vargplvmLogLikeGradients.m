function g = vargplvmLogLikeGradients(model)

% VARGPLVMLOGLIKEGRADIENTS Compute the gradients for the variational GPLVM.
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
% COPYRIGHT : Michalis K. Titsias, 2009, 2010,2011
%
% COPYRIGHT :  Mauricio Alvarez, 2009, 2010,2011
%
% COPYRIGHT :  Neil D. Lawrence, 2009, 2010,2011
%
% COPYRIGHT : Andreas Damianou, 2010,2011

% SEEALSO : vargplvmLogLikelihood, vargplvmCreate, modelLogLikeGradients

% VARGPLVM


% The gradient of the kernel of the dynamics (e.g. temporal prior)
gDynKern = [];

% KL divergence terms: If there are no dynamics, the only params w.r.t to
% which we need derivatives are the var. means and covars. With dynamics
% (see below) we also need derivs. w.r.t theta_t (dyn. kernel's params).
if ~isfield(model, 'dynamics') || isempty(model.dynamics)
    if ~(isfield(model, 'onlyLikelihood') && model.onlyLikelihood)
        gVarmeansKL = - model.vardist.means(:)';
        % !!! the covars are optimized in the log space (otherwise the *
        % becomes a / and the signs change)
        gVarcovsKL = 0.5 - 0.5*model.vardist.covars(:)';
    else
        gVarmeansKL = 0;
        gVarcovsKL = 0;
    end
end

%fprintf(' %d \n',sum([gVarmeansKL gVarcovsKL]));%%%%%TEMP

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


% if strcmp(model.kern.type, 'rbfardjit')
%     % different derivatives for the variance, which is super-numerically stable for
%     % this particular kernel
%     gKern(1) = 0.5*model.d*( - model.k+ sum(sum(model.invLat.*model.invLat))/model.beta - model.beta*(model.Psi0-model.TrC)  )...
%                     + 0.5*tmpV;
% end

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
    
    % This should work much faster
    %gX_u2 = kernKuuXuGradient(model.kern, model.X_u, gK_uu);
    
    %sum(abs(gX_u2(:)-gX_u(:)))
    %pause
    
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


% If we only want to exclude the derivatives for the variational
% distribution, the following big block will be skipped.
if ~(isfield(model, 'onlyKernel') && model.onlyKernel)
    if isfield(model, 'dynamics') && ~isempty(model.dynamics)
        % Calculate the derivatives for the reparametrized variational and Kt parameters.
        % The formulae for these include in a mixed way the derivatives of the KL
        % term w.r.t these., so gVarmeansKL and gVarcovsKL are not needed now. Also
        % the derivatives w.r.t kernel parameters also require the derivatives of
        % the likelihood term w.r.t the var. parameters, so this call must be put
        % in this part.
        
        % For the dynamical GPLVM further the original covs. must be fed,
        % before amending with the partial derivative due to exponing to enforce
        % positiveness.
        gVarcovsLik = gVarcovs0 + gVarcovs1 + gVarcovs2;
        [gVarmeans gVarcovs gDynKern] = modelPriorReparamGrads(model.dynamics, gVarmeansLik, gVarcovsLik, gIndRep);
        % Variational variances are positive: Now that the final covariances
        % are obtained we amend with the partial derivatives due to the
        % exponential transformation to ensure positiveness.
         if ~isfield(model, 'notransform') || (isfield(model,'notransform') && model.notransform == false)
       		 gVarcovs = (gVarcovs(:).*model.dynamics.vardist.covars(:))';
    	end
    else
        % For the non-dynamical GPLVM these cov. derivatives are the final, so
        % it is time to amend with the partial derivative due to exponing them
        % to force posigiveness.
        gVarcovs0 = (gVarcovs0(:).*model.vardist.covars(:))';
        gVarcovs1 = (gVarcovs1(:).*model.vardist.covars(:))';
        gVarcovs2 = (gVarcovs2(:).*model.vardist.covars(:))';
        
        gVarcovsLik = gVarcovs0 + gVarcovs1 + gVarcovs2;
        gVarmeans = gVarmeansLik + gVarmeansKL;
        %gVarcovsLik = (gVarcovsLik(:).*model.vardist.covars(:))';
        gVarcovs = gVarcovsLik + gVarcovsKL;
    end
else
    gVarmeans = [];
    gVarcovs = [];
    gDynKern = [];
end



%{
    % %%% Compute Gradients with respect to X_u %%%
    % gKX = kernGradX(model.kern, model.X_u, model.X_u);
    %
    % % The 2 accounts for the fact that covGrad is symmetric
    % gKX = gKX*2;
    % dgKX = kernDiagGradX(model.kern, model.X_u);
    % for i = 1:model.k
    %   gKX(i, :, i) = dgKX(i, :);
    % end
    %
    % % Allocate space for gX_u
    % gX_u = zeros(model.k, model.q);
    % % Compute portion associated with gK_u
    % for i = 1:model.k
    %   for j = 1:model.q
    %     gX_u(i, j) = gKX(:, j, i)'*gK_uu(:, i);
    %   end
    % end
    %
    % % This should work much faster
    % %gX_u2 = kernKuuXuGradient(model.kern, model.X_u, gK_uu);
    %
    % %sum(abs(gX_u2(:)-gX_u(:)))
    % %pause
    %
    % gInd = gInd1 + gInd2 + gX_u(:)';
%}



%%% TEMP:  the following needs some more testing...
% If fixInducing is true then the inducing points are not optimised but are
% rather reparametrized as X_u=X=mu (static case). Then, there is no derivative
% for X_u any more but the one for mu must be amended by adding the partial
% derivatives of the inducing points (which is just gInd because X_u is
% mapped to mu with the identity function whose deriv. is 1). In the
% dynamics case X_u=X=mu=Kt*mu_bar so we further need to amend with
% d mu/ d mu_bar = K_t because we need
% d F/ d mu_bar instead of d F/ d mu.
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

% gVarmeans = gVarmeans0 + gVarmeans1 + gVarmeans2 + gVarmeansKL;
% gVarcovs = gVarcovs0 + gVarcovs1 + gVarcovs2 + gVarcovsKL;



if isfield(model.vardist,'paramGroups')
    gVar = gVar*model.vardist.paramGroups;
end


% If we only want to exclude the derivatives for the variational
% distribution, the following big block will be skipped.
if ~(isfield(model, 'onlyKernel') && model.onlyKernel)
    % It may better to start optimize beta a bit later so that
    % so that the rest parameters can be initialized
    % (this could help to escape from the trivial local
    % minima where the noise beta explains all the data)
    % The learnBeta option deals with the above.
    
    % This constrains the variance of the dynamics kernel to one
    % (This piece of code needs to be done in better way with unit variance dynamic
    %  kernels. The code below also will only work for rbf dynamic kernel)
    % Assume that the base rbf/matern/etc kernel is first in the compound
    % structure
    if isfield(model, 'dynamics') && ~isempty(model.dynamics)
        if isfield(model.dynamics.kern, 'comp')
            if strcmp(model.dynamics.kern.comp{1}.type,'rbf') || ...
                    strcmp(model.dynamics.kern.comp{1}.type,'matern32') || ...
                    strcmp(model.dynamics.kern.comp{1}.type,'rbfperiodic') || ...
                    strcmp(model.dynamics.kern.comp{1}.type,'rbfperiodic2')
                if ~isfield(model.dynamics, 'learnVariance') || ~model.dynamics.learnVariance
                    gDynKern(2) = 0;
                end
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
        
        %--
       % if strcmp(model.dynamics.kern.type,'invcmpnd')
            if ~isfield(model.dynamics, 'learnVariance') || ~model.dynamics.learnVariance
                % The field model.dynamics.fixedVariance must have the
                % indexes of the gradient vector that correspond to
                % variances of kernels of the matern class and we wish to
                % not learn them.
                if isfield(model.dynamics, 'fixedKernVariance') && ~isempty(model.dynamics.fixedKernVariance)
                    gDynKern(model.dynamics.fixedKernVariance) = 0;
                end
            end
        % end
        %--
    end
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

if ~learnInducing
    gInd = [];
end

% At this point, gDynKern will be [] if there are no dynamics.
g = [gVar gDynKern gInd gKern gBetaFinal];


% A 'dirty' trick to fix some parameters
if isfield(model, 'fixParamIndices') && ~isempty(model.fixParamIndices)
    g(model.fixParamIndices) = 0;
end

% if model.learnBeta == 1
%     % optimize all parameters including beta
%     g = [gVar gDynKern gInd gKern gBeta];
% else
%     % keep beta fixed
%     g = [gVar gDynKern gInd gKern 0*gBeta];
% end


%if model.kern.comp{1}.variance > 100 | model.kern.comp{1}.variance < 0.001
%     model.kern.comp{1}
%end

% delete afterwards
%g = [gVarmeans2 gVarcovs2 (gX_u(:)'+ gInd2) (gKern3+gKern2) 0*gBeta];

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

%%%%TEMP
%{
    load TEMPbetaGradTrC;
    TEMPbetaGradTrC = [TEMPbetaGradTrC model.d*0.5*model.TrC];
    save 'TEMPbetaGradTrC.mat' 'TEMPbetaGradTrC';

    load TEMPbetaGradNksigm;
    TEMPbetaGradNksigm=[TEMPbetaGradNksigm model.d*0.5*(model.N-model.k)*sigm];
    save 'TEMPbetaGradNksigm.mat' 'TEMPbetaGradNksigm';

    load TEMPbetaGradPsi0;
    TEMPbetaGradPsi0=[TEMPbetaGradPsi0 (-0.5*model.d*model.Psi0)];
    save 'TEMPbetaGradPsi0.mat' 'TEMPbetaGradPsi0';

    load TEMPbetaGradTrPP;
    TEMPbetaGradTrPP=[TEMPbetaGradTrPP 0.5*model.TrPP];
    save 'TEMPbetaGradTrPP.mat' 'TEMPbetaGradTrPP';

    load TEMPbetaGradLat;
    TEMPbetaGradLat=[TEMPbetaGradLat (1/(model.beta^2)) * model.d * sum(sum(model.invLat.*model.invLat))*0.5];
    save 'TEMPbetaGradLat.mat' 'TEMPbetaGradLat';

    load TEMPbetaGradPlm;
    TEMPbetaGradPlm=[TEMPbetaGradPlm sigm*sum(sum(PLm.*PLm))*0.5];
    save 'TEMPbetaGradPlm.mat' 'TEMPbetaGradPlm';
%}
%%%%%


%gBeta = 0.5*(model.d*(model.TrC + (model.N-model.k)*sigm -model.Psi0) ...
%	- model.TrYY + model.TrPP ...
%	+ sigm * sum(sum(model.K_uu .* model.Tb)));

if ~isstruct(model.betaTransform)
    fhandle = str2func([model.betaTransform 'Transform']);
    gBeta = gBeta*fhandle(model.beta, 'gradfact');
else
    fhandle = str2func([model.betaTransform.type 'Transform']);
    gBeta = gBeta*fhandle(model.beta, 'gradfact', model.betaTransform.transformsettings);
end


g_Lambda = repmat(-0.5*model.beta*model.d, 1, model.N);
end
