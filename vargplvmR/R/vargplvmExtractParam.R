vargplvmExtractParam <-
function (model, only.values = TRUE, untransformed.values = FALSE)
{
# % VARGPLVMEXTRACTPARAM Extract a parameter vector from a variational GP-LVM model.
# % FORMAT
# % DESC extracts a parameter vector from a given VARGPLVM structure.
# % ARG model : the model from which parameters are to be extracted.
# % RETURN params : the parameter vector extracted from the model.
# %
# % DESC does the same as above, but also returns parameter names.
# % ARG model2 : the model structure containing the information about
# % the model.
# % RETURN params : a vector of parameters from the model.
# % RETURN names : cell array of parameter names.
# %
# % COPYRIGHT : Michalis K. Titsias and Neil D. Lawrence, 2009-2011
# %
# % Modifications: Andreas C. Damianou, 2010-2011 
# %
# % SEEALSO : vargplvmCreate, vargplvmExpandParam, modelExtractParam
# 
# % VARGPLVM

# % Parameters must be returned as a vector in the following order (left to right) 
# % - parameter{size} -
# % vardistParams{model.vardist.nParams} % mu, S
# %       OR
# % [dynamicsVardistParams{dynamics.vardist.nParams} dynamics.kernParams{dynamics.kern.nParams}] % mu_bar, lambda
# % inducingInputs{model.q*model.k}
# % kernelParams{model.kern.nParams}
# % beta{prod(size(model.beta))}


if (!only.values)
returnNames <- TRUE 
else
returnNames <- FALSE 


if (("dynamics" %in% names(model)) && !is.null(model$dynamics))
{
  cat("vargplvmExtractParam to do ..\n")
# % [VariationalParameters(reparam)   dynKernelParameters]
  #   To do
  # # if (returnNames)
  # [dynParams, dynParamNames] <- modelExtractParam(model.dynamics) 
  # names <- dynParamNames 
  # else
  #   dynParams <- modelExtractParam(model.dynamics) 
  # end
  # params <- dynParams 
} else {
# % Variational parameters 
  if (returnNames)
  {
    cat("to do vargplvmExtractParam\n")
# %[varParams, varNames] <- vardistExtractParam(model.vardist) 
    # [varParams, varNames] <- modelExtractParam(model$vardist) 
# %names <- varNames{:}  %%% ORIGINAL
    # names <- varNames  %%% NEW
  }
  else {
# %varParams <- vardistExtractParam(model$vardist) 
    varParams <- modelExtractParam(model$vardist) 
  }
params <- varParams 
}


# % Inducing inputs 
if (!model$fixInducing)
{
  params <- c(params, c(matrix(model$X_u, 1, prod(dim(model$X_u)))))
  if (returnNames) 
  {
    for (i in 1:dim(model$X_u)[1])
    {
      for (j in 1:dim(model$X_u)[2])
        X_uNames[i, j] <- paste("X_u(", i, ",", j, ")", sep = "") 
    }
    names <- c(matrix(names, ncol = 1) , matrix(X_uNames, ncol = 1)) 
  } 
}
#     % Kernel parameters  
#     if returnNames
#     [kernParams, kernParamNames] <- kernExtractParam(model$kern)  
#     for i <- 1:length(kernParamNames)
#     kernParamNames{i} <- ['Kernel, ' kernParamNames{i}] 
#     end
#     names <- {names{:}, kernParamNames{:}} 
#     else
      kernParams <- kernExtractParam(model$kern, matlabway = TRUE) 
#     end
    params <- c(params, kernParams) 
    
    
#     % beta in the likelihood 
if (model$optimiseBeta)
{
    if (!is.list(model$betaTransform))
    {
      fhandle <- paste(model$betaTransform, "Transform", sep = "") 
      betaParam <- do.call(fhandle, list(model$beta, "xtoa", matlabway = TRUE)) 
    } else {
      if (("transformsettings" %in% names(model$betaTransform)) && (length(model$betaTransform$transformsettings)>0))
      {
        fhandle <- paste(model$betaTransform$type, "Transform", sep = "") 
        betaParam <- do.call(fhandle, list(model$beta, "xtoa", model$betaTransform$transformsettings)) 
      } else {
        stop("vargplvmExtractParam: Invalid transform specified for beta.\n")
      }
    }    
    params <- c(params, matrix(betaParam, nrow = 1))
   
    if (returnNames)
    {
      betaParamNames <- paste("Beta", 1:length(betaParam), sep="") 
      names(params) <- c(names, betaParamNames)
    }
}
return (params)
}
