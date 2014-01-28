vargplvmPosteriorVar <-
function (model, X)
{
# % VARGPLVMPOSTERIORVAR variances of the posterior at points given by X.
# % FORMAT
# % DESC returns the posterior mean and variance for a given set of
# % points.
# % ARG model : the model for which the posterior will be computed.
# % ARG x : the input positions for which the posterior will be
# % computed.
# % RETURN mu : the mean of the posterior distribution.
# % RETURN sigma : the variances of the posterior distributions.
# %
# % SEEALSO : gpPosteriorMeanVar, vargplvmCreate
# %
# % COPYRIGHT : Neil D. Lawrence, 2005, 2006
# 
# % VARGPLVM


# %% ORIGINAL
model$K_uf <- kernCompute(model$kern, model$X_u, X) 
model$A <- (1/model$beta)*model$K_uu + model$K_uf%*%t(model$K_uf) 
temp <- pdinv(model$A) 
model$Ainv <- temp$Ainv
U <- temp$UC
varsigma <- gpPosteriorVar(model, X) 

# %%% TEMP
# %[void, varsigma] <- vargplvmPosteriorMeanVar(model, X) 
return (varsigma)
}
