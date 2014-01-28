vargplvmLogLikeGradients <-
function (model)
{
# % VARGPLVMLOGLIKEGRADIENTS Compute the gradients for the variational GPLVM.
# % FORMAT
# % DESC returns the gradients of the log likelihood with respect to the
# % parameters of the GP-LVM model and with respect to the latent
# % positions of the GP-LVM model.
# % ARG model : the FGPLVM structure containing the parameters and
# % the latent positions.
# % RETURN g : the gradients of the latent positions (or the back
# % constraint's parameters) and the parameters of the GP-LVM model.
# %
# % FORMAT
# % DESC returns the gradients of the log likelihood with respect to the
# % parameters of the GP-LVM model and with respect to the latent
# % positions of the GP-LVM model in seperate matrices.
# % ARG model : the FGPLVM structure containing the parameters and
# % the latent positions.
# % RETURN gX : the gradients of the latent positions (or the back
# % constraint's parameters).
# % RETURN gParam : gradients of the parameters of the GP-LVM model.
# %
# % COPYRIGHT : Michalis K. Titsias, 2009, 2010,2011
# %
# % COPYRIGHT :  Mauricio Alvarez, 2009, 2010,2011
# %
# % COPYRIGHT :  Neil D. Lawrence, 2009, 2010,2011
# %
# % COPYRIGHT : Andreas Damianou, 2010,2011
# 
# % SEEALSO : vargplvmLogLikelihood, vargplvmCreate, modelLogLikeGradients
# 
# % VARGPLVM


# % The gradient of the kernel of the dynamics (e.g. temporal prior)
gDynKern <- NULL

# % KL divergence terms: If there are no dynamics, the only params w.r.t to
# % which we need derivatives are the var. means and covars. With dynamics
# % (see below) we also need derivs. w.r.t theta_t (dyn. kernel's params).
if ((!("dynamics" %in% names(model))) || is.null(model$dynamics))
{
    if (!((("onlyLikelihood" %in% names(model))) && model$onlyLikelihood))
    {
        gVarmeansKL <- - matrix(model$vardist$means, nrow = 1)
#         % !!! the covars are optimized in the log space
        gVarcovsKL <- 0.5 - 0.5*matrix(model$vardist$covars, nrow = 1) 
    } else {
        gVarmeansKL <- 0 
        gVarcovsKL <- 0 
    }
}

# %fprintf(' %d \n',sum([gVarmeansKL gVarcovsKL])) %%%%%TEMP

# % Likelihood terms (coefficients)
likhod <- vargpCovGrads(model) 
# [gK_uu, gPsi0, gPsi1, gPsi2, g_Lambda, gBeta, tmpV]

# % Get (in three steps because the formula has three terms) the gradients of
# % the likelihood part w.r.t the data kernel parameters, variational means
# % and covariances (original ones). From the field model.vardist, only
# % vardist.means and vardist.covars and vardist.lantentDimension are used.

kvp1g <- kernVardistPsi1Gradient(model$kern, model$vardist, model$X_u, t(likhod$gPsi1)) 

# [gKern1, gVarmeans1, gVarcovs1, gInd1]                                                
kvp2g <- kernVardistPsi2Gradient(model$kern, model$vardist, model$X_u, likhod$gPsi2) 
# [gKern2, gVarmeans2, gVarcovs2, gInd2] 

kvp0g <- kernVardistPsi0Gradient(model$kern, model$vardist, likhod$gPsi0) 
# [gKern0, gVarmeans0, gVarcovs0]             

gKern3 <- kernGradient(model$kern, model$X_u, matlabway = TRUE, likhod$gK_uu) 

# % At this point, gKern gVarmeansLik and gVarcovsLik have the derivatives for the
# % likelihood part. Sum all of them to obtain the final result.
gKern <- kvp0g$gKern + kvp1g$gKern + kvp2g$gKern + gKern3 
gVarmeansLik <- kvp0g$gVarmeans + kvp1g$gVarmeans + kvp2g$gVarmeans 

if (model$kern$type == "rbfardjit")
{
#     % different derivatives for the variance, which is super-numerically stable for
#     % this particular kernel
  cat("to do vLoglikeG")
#     if (model$learnSigmaf == 1)
#     {
#         gKern[1] <- 0.5*model$d*( - model$k+ sum(sum(model$invLat.*model$invLat))/model$beta - model$beta*(model$Psi0-model$TrC)  )...
#             + 0.5*likhod$tmpV 
#         
#         if ~isstruct(model$kern$transforms(1))
#             fhandle <- str2func([model$kern$transform(1) 'Transform']) 
#             gKern(1) <- gKern(1).*fhandle(model$kern$variance, 'gradfact') 
#         else
#             fhandle <- str2func([model$kern$transforms(1)$type 'Transform']) 
#             if ~isfield(model$kern$transforms(1), 'transformsettings')
#                 gKern(1) <- gKern(1).*fhandle(model$kern$variance, 'gradfact') 
#             else
#                 gKern(1) <- gKern(1).*fhandle(model$kern$variance, 'gradfact', model$kern$transforms(1)$transformsettings) 
#             end
#         end
#     } else {
#         gKern(1) <- 0 
#     }
}

# %%% Compute Gradients with respect to X_u %%%
gKX <- kernGradX(model$kern, model$X_u, model$X_u, matlabway = TRUE) 

# % The 2 accounts for the fact that covGrad is symmetric
gKX <- gKX*2 

dgKX <- kernDiagGradX(model$kern, model$X_u) 

for (i in 1:model$k)
    gKX[i, , i] <- dgKX[i, ] 


# % Allocate space for gX_u
gX_u <- matrix(0, model$k, model$q) 
# % Compute portion associated with gK_u

for (i in 1:model$k)
{
    for (j in 1:model$q)
        gX_u[i, j] <- t(gKX[, j, i])%*%likhod$gK_uu[, i] 
}
# % This should work much faster
# %gX_u2 <- kernKuuXuGradient(model$kern, model$X_u, gK_uu) 
# 
# %sum(abs(gX_u2(:)-gX_u(:)))
# %pause

gInd <- kvp1g$gInd + kvp2g$gInd + matrix(gX_u, nrow = 1) 

# % If the inducing points are fixed (tied to the latent points) then
# % X_u<-K_t*dynamics$vardist$means and the derivatives w.r.t theta_t must be
# % amended with the appropriate partial derivatives. gInd must be passed,
# % in that case, as an argument to the function which calculates the
# % derivatives for the reparametrized quantities.
if (("fixInducing" %in% names(model)) && model$fixInducing)
    gIndRep <- gInd 
else
    gIndRep<-NULL 

# % If we only want to exclude the derivatives for the variational
# % distribution, the following big block will be skipped.
if (!(("onlyKernel" %in% names(model)) && model$onlyKernel))
{
    if (("dynamics" %in% names(model)) && !is.null(model$dynamics))
    {
      cat("to do vpLLG")
#         % Calculate the derivatives for the reparametrized variational and Kt parameters.
#         % The formulae for these include in a mixed way the derivatives of the KL
#         % term w.r.t these., so gVarmeansKL and gVarcovsKL are not needed now. Also
#         % the derivatives w.r.t kernel parameters also require the derivatives of
#         % the likelihood term w.r.t the var. parameters, so this call must be put
#         % in this part.
#         
#         % For the dynamical GPLVM further the original covs. must be fed,
#         % before amending with the partial derivative due to exponing to enforce
#         % positiveness.
#         gVarcovsLik <- gVarcovs0 + gVarcovs1 + gVarcovs2 
#         [gVarmeans gVarcovs gDynKern] <- modelPriorReparamGrads(model$dynamics, gVarmeansLik, gVarcovsLik, gIndRep) 
#         % Variational variances are positive: Now that the final covariances
#         % are obtained we amend with the partial derivatives due to the
#         % exponential transformation to ensure positiveness.
#          if ~isfield(model, 'notransform') || (isfield(model,'notransform') && model$notransform == false)
#           gVarcovs <- (gVarcovs(:).*model$dynamics$vardist$covars(:))' 
#      end
    } else {
#         % For the non-dynamical GPLVM these cov. derivatives are the final, so
#         % it is time to amend with the partial derivative due to exponing them
#         % to force posigiveness.
        gVarcovs0 <- (matrix(kvp0g$gVarcovars, nrow = 1)*matrix(model$vardist$covars, nrow = 1)) 
        gVarcovs1 <- (matrix(kvp1g$gVarcovars, nrow = 1)*matrix(model$vardist$covars, nrow = 1)) 
        gVarcovs2 <- (matrix(kvp2g$gVarcovars, nrow = 1)*matrix(model$vardist$covars, nrow = 1)) 
        
        gVarcovsLik <- gVarcovs0 + gVarcovs1 + gVarcovs2 
        gVarmeans <- gVarmeansLik + gVarmeansKL 
#         %gVarcovsLik <- (gVarcovsLik(:).*model$vardist$covars(:))' 
        gVarcovs <- gVarcovsLik + gVarcovsKL 
    }
} else {
    gVarmeans <- NULL 
    gVarcovs <- NULL
    gDynKern <- NULL
}

# %%% TEMP:  the following needs some more testing...
# % If fixInducing is true then the inducing points are not optimised but are
# % rather reparametrized as X_u=X=mu (static case). Then, there is no derivative
# % for X_u any more but the one for mu must be amended by adding the partial
# % derivatives of the inducing points (which is just gInd because X_u is
# % mapped to mu with the identity function whose deriv. is 1). In the
# % dynamics case X_u=X=mu=Kt*mu_bar so we further need to amend with
# % d mu/ d mu_bar = K_t because we need
# % d F/ d mu_bar instead of d F/ d mu.
if (("fixInducing" %in% names(model)) && model$fixInducing)
{
# % If there are dynamics the derivative must further be amended with the
#     % partial deriv. due to the mean reparametrization.
  cat("to do vgplvmLLG fixInducing")
#     if isfield(model, 'dynamics') && ~isempty(model$dynamics)
#         gInd <- reshape(gInd,model$k,model$q) 
#         %gInd <- gInd' * model$dynamics$Kt 
#         gInd <-  model$dynamics$Kt * gInd 
#         gInd <- gInd(:)' 
#     end
#     %gVarmeans(model$inducingIndices, :) <- gVarmeans(model$inducingIndices,
#     %:) + gInd  % This should work AFTER reshaping the matrices...but here
#     %we use all the indices anyway.
#     gVarmeans <- gVarmeans + gInd 
#     gInd <- []  % Inducing points are not free variables anymore, they dont have derivatives on their own.
}

gVar <- c(gVarmeans, gVarcovs) 

# % gVarmeans <- gVarmeans0 + gVarmeans1 + gVarmeans2 + gVarmeansKL 
# % gVarcovs <- gVarcovs0 + gVarcovs1 + gVarcovs2 + gVarcovsKL 

if ("paramGroups" %in% names(model$vardist))
    gVar <- gVar*model$vardist$paramGroups 


# % If we only want to exclude the derivatives for the variational
# % distribution, the following big block will be skipped.
if (!(("onlyKernel" %in% names(model)) && model$onlyKernel))
{
#     % It may better to start optimize beta a bit later so that
#     % so that the rest parameters can be initialized
#     % (this could help to escape from the trivial local
#     % minima where the noise beta explains all the data)
#     % The learnBeta option deals with the above.
#     
#     % This constrains the variance of the dynamics kernel to one
#     % (This piece of code needs to be done in better way with unit variance dynamic
#     %  kernels. The code below also will only work for rbf dynamic kernel)
#     % Assume that the base rbf/matern/etc kernel is first in the compound
#     % structure
    if (("dynamics" %in% names(model)) && !is.null(model$dynamics))
    {
      cat("to do vpglvmLLG")
#         if strcmp(model$dynamics$kern$comp{1}$type,'rbf') || strcmp(model$dynamics$kern$comp{1}$type,'matern32') || strcmp(model$dynamics$kern$comp{1}$type,'rbfperiodic') || strcmp(model$dynamics$kern$comp{1}$type,'rbfperiodic2')
#             if ~isfield(model$dynamics, 'learnVariance') || ~model$dynamics$learnVariance
#                 gDynKern(2) <- 0 
#             end
#         end
#         
#         
#         %___NEW: assume that the second rbf/matern etc kernel is last in the
#         %compound kernel
#         %if numel(model$dynamics$kern$comp) > 3
#         if isfield(model$dynamics, 'learnSecondVariance') && ~model$dynamics$learnSecondVariance   %%%%% NEW
#             gDynKern(end) <- 0 
#         end
#         %end
#         %___
    }
}

# % In case we are in the phase where the vardistr. is initialised (see above
# % for the variance of the kernel), beta is kept fixed. For backwards
# % compatibility this can be controlled either with the learnBeta field or
# % with the initVardist field. The later overrides the first.
if (("learnBeta" %in% names(model)) && model$learnBeta)
    gBetaFinal <- likhod$gBeta 
else
    gBetaFinal <- 0*likhod$gBeta 

if ("initVardist" %in% names(model))
{
    if (model$initVardist == 1)
        gBetaFinal <- 0*likhod$gBeta 
    else
        gBetaFinal <- likhod$gBeta 
}

# % At this point, gDynKern will be [] if there are no dynamics.
g <- c(gVar, gDynKern, gInd, gKern, gBetaFinal) 

return (g)
}
