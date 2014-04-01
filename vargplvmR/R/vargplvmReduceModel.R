vargplvmReduceModel <-
function (model, P, dims, only.values =  TRUE)
{
# % VARGPLVMREDUCEMODEL prunes out dimensions of the model.
# % FORMAT
# % DESC order the latent dimensions acorrding to the inputScales and
# % reduces the model to have smaller number of latent dimensions.
# % ARG model : the model to be reduced.
# % ARG P : the number of dimensions to move to (setting to model.q will
# % just reorder the dimensions in the model).
# % ARG dims : (optional) explicit set of dimensions to use
# % RETURN model : the model with the reduced number of dimensions.
# %
# % COPYRIGHT : Michalis K. Titsias, 2009
# %  
# % MODIFICATIONS : Neil D. Lawrence, 2009, Patrick Sauer, 2011
# % 
# % SEEALSO : vargplvmCreate  
  # 
# % VARGPLVM
  
  if (nargs() == 3)
    P <- length(dims)
  
  
# % create temporary model
  options <- vargplvmOptions("dtcvar")
  options$kern <- NULL
  options$numActive <- model$k 
  if (!(model$kern$type == "cmpnd"))
  {
    options$kern <- model$kern$type 
  } else if (model$kern$type  == "cmpnd") {
    options$kern[1] <- model$kern$comp[[1]]$type 
    for (i in 2:length(model$kern$comp))
      options$kern[i] <- model$kern$comp[[i]]$type 
  }
  
  if (is.list(model$betaTransform))
    options$betaTransform <- model$betaTransform  
  
  mm <- vargplvmCreate(P, model$d, model$y, options) 
  N <- dim(model$vardist$means)[1] 
  
  if (!(model$kern$type == "cmpnd"))
  { 
    #   to do 
    cat("do to vargplvmReduceModel \n")
    #   if ((model$kern$type == "rbfardjit") || (model$kern$type == "linard2") || (model$kern$type == "rbfard2"))
    #   {
    #     if (nargs() == 2)
    #     {
    #         order2 <- order(-model$kern$inputScales) 
    #         mm$kern$inputScales <- model$kern$inputScales(order2[1:P]) 
    #     } else {
    #         order2 <- c(dims, setdiff(1:length(model$kern$inputScales), dims)] 
    #         mm$kern$inputScales <- model$kern$inputScales(order2) 
    #     }   
    #   }
  }
  else
  {
    for (i in 1:length(model$kern$comp))
    {
      if ((model$kern$comp[[i]]$type == "rbfardjit") || 
        (model$kern$comp[[i]]$type == "linard2") || 
        (model$kern$comp[[i]]$type == "rbfard2"))
      {  
        if (nargs() == 2)
        {
          order2 <- order(-model$kern$comp[[i]]$inputScales) 
          mm$kern$comp[[i]]$inputScales <- model$kern$comp[[i]]$inputScales[order2[1:P]] 
        } else {
          order2 <- c(dims, setdiff(1:length(model$kern$comp[[i]]$inputScales), dims)) 
          mm$kern$comp[[i]]$inputScales <- model$kern$comp[[i]]$inputScales[order2] 
        }
        #       % you order only wrt the first ARD kernel you find 
        break   
      }
    }
  }
  
  mm$vardist$means  <- model$vardist$means[,order2[1:P]] 
  mm$vardist$covars <- model$vardist$covars[,order2[1:P]] 
  
  mm$X_u <- model$X_u[,order2[1:P]]
  mm$X   <- model$vardist$means[,order2[1:P]] 
  mm$inputSclsOrder <- order2 
  
  initParams   <- vargplvmExtractParam(mm) 
  mm$numParams <- length(initParams) 
  
# % This forces kernel computation.
  mm <- vargplvmExpandParam(mm, initParams) 
  if (only.values)
    return (mm)
  else
    return (list(mm = mm, order = order2))
}
