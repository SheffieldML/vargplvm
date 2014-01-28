vargplvmParamInit <-
  function (model,Y,X)
  {
# % VARGPLVMPARAMINIT Initialize the variational GPLVM from the data
# % FORMAT
# % DESC no description.
# % COPYRIGHT: Michalis Titsias 2009-2011
# % VARGPLVM
    
#       Y = model$m 
#       X = model$X
    if (!(model$kern$type == "cmpnd"))
    {   
      if ((model$kern$type == "rbfard2")  || (model$kern$type == "rbfardjit"))
      {
        model$kern$inputScales <- 5/((apply(X, 2, max) -apply(X, 2, min))^2) 
        #        %model$kern$variance <- max(var(Y))  % !!! better use mean(var(Y)) here...
        model$kern$variance <- mean(apply(Y, 2, var))  #% NEW!!!
      }  else if (model$kern$type == "linard2") 
      {
        model$kern$inputScales <- 5/((apply(X, 2, max) -apply(X, 2, min))^2) 
      }
      
    } else {
      
      for (i in 1:length(model$kern$comp))
      {
        if ((model$kern$comp[[i]]$type == "rbfard2") || (model$kern$type == "rbfardjit")) 
        {
          model$kern$comp[[i]]$inputScales <- 5/((apply(X, 2, max) -apply(X, 2, min))^2) 
          #          %model$kern$comp{i}$variance <- max(var(Y)) 
          model$kern$comp[[i]]$variance <- var(matrix(model$m,dim(model$m)[1]*dim(model$m)[2],1)) 
        } else if (model$kern$comp[[i]]$type == "linard2")
        {
          model$kern$comp[[i]]$inputScales <- 0.01*max(apply(Y, 2, var))*matrix(1, 1,dim(X)[2]) #% %5./(((max(X)-min(X))).^2)   
        }  
      }
    }
     
    model$beta <- 1000 #%/max(var(Y)) 
    
    
    initParams <- vargplvmExtractParam(model) 
    model$numParams <- length(initParams) 
# % This forces kernel computation.
    model <- vargplvmExpandParam(model, initParams) 
    return (model)
  }
