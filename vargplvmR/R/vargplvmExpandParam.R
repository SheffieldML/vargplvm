vargplvmExpandParam <-
function (model, params)
{
# % VARGPLVMEXPANDPARAM Expand a parameter vector into a GP-LVM model.
# % FORMAT
# % DESC takes an VARGPLVM structure and a vector of parameters, and
# % fills the structure with the given parameters. Also performs any
# % necessary precomputation for likelihood and gradient
# % computations, so can be computationally intensive to call.
# % ARG model : the VARGPLVM structure to put the parameters in.
# % ARG params : parameter vector containing the parameters to put in
# % the VARGPLVM structure.
# % 
# %
# % COPYRIGHT : Michalis K. Titsias, 2009-2011
# % COPYRIGHT : Neil D. Lawrence, 2009-2011
# % 
# % Modifications: Andreas C. Damianou, 2010-2011 
# %
# % SEEALSO : vargplvmCreate, vargplvmExtractParam, modelExpandParam
# 
# % VARGPLVM

startVal <- 1 
if (("dynamics" %in% names(model)) && (length(model$dynamics) > 0))
{
#     % variational parameters (reparametrized) AND dyn$kernel's parameters
    endVal <- model$dynamics$nParams 
    model$dynamics <- modelExpandParam(model$dynamics, params(startVal:endVal)) 
} else {
#     % variational parameters (means and covariances), original ones
    endVal <- model$vardist$nParams 
    model$vardist <- modelExpandParam(model$vardist, params[startVal:endVal])  
}

# % inducing inputs 
startVal <- endVal+1 
if (model$fixInducing)
{
  cat("to do vargplvmExpandParam\n")
  # TO DO
#     if (("dynamics" %in% names(model)) && (length(model$dynamics) > 0))
#     {
#         Kt <- kernCompute(model$dynamics$kern, model$dynamics$t) #%%%%%%%%%%%%5
#         model$X_u <- Kt*model$dynamics$vardist$means  %dynamics
#         model$X_u <- model$X_u(model$inducingIndices,:) 
#     } else {
#         model$X_u <- model$vardist$means(model$inducingIndices, :)  % static
#     }
#     % X_u values are taken from X values.
#    % model$X_u <- model$X(model$inducingIndices, :) 
} else {
#     % Parameters include inducing variables.
    endVal <- endVal + model$q*model$k 
    model$X_u <- matrix(params[startVal:endVal], model$k, model$q) 
}

#     % kernel hyperparameters 
startVal <- endVal+1  
endVal <- endVal + model$kern$nParams 
model$kern <- kernExpandParam(model$kern, params[startVal:endVal], 
matlabway = TRUE) 

# % likelihood beta parameters
if (model$optimiseBeta)
{
  startVal <- endVal + 1 
  endVal <- endVal + prod(dim(model$beta)) 
  if (!is.list(model$betaTransform))
  {
    fhandle <- paste(model$betaTransform, "Transform", sep ="")
    model$beta <- do.call(fhandle, list(params[startVal:endVal], "atox",
     matlabway = TRUE))
  } else {
      if (("transformsettings" %in% names(model$betaTransform)) && (length(model$betaTransform$transformsettings) > 0))
      {
          fhandle <- paste(model$betaTransform$type, "Transform", sep = "")
          model$beta <- do.call(fhandle, list(params[startVal:endVal], "atox", model$betaTransform$transformsettings))
      } else {
          stop ("vargplvmExtractParam: Invalid transform specified for beta.") 
      }
  }
}

model$nParams <- endVal 

# % Update statistics
model <- vargplvmUpdateStats(model, model$X_u) 

# % %%%TEMP: This is not needed, probably. If yes, it should be merged with
# % %%%the above code for fixInducing.
# % if model$fixInducing
# %     model$X_u<-model$X(model$inducingIndices, :) 
# % end
# % %%%
return (model)
}
