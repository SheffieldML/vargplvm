gpPosteriorVar <-
  function (model, X) 
  {
# % GPPOSTERIORVAR Variances of the posterior at points given by X.
# % FORMAT
# % DESC returns the posterior variance for a given set of
# % points.
# % ARG model : the model for which the posterior will be computed.
# % ARG x : the input positions for which the posterior will be
# % computed.
# % RETURN sigma : the variances of the posterior distributions.
# %
# % SEEALSO : gpCreate, gpPosteriorMeanCovar, gpPosteriorGradMeanVar, gpPosteriorMeanVar
# %
# % COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2009
    # 
# % GP
    
    
    maxMemory <-  1000000 
    switch (EXPR = model$approx,
            ftc = {
              chunkSize <-  ceiling(maxMemory/model$N)},
            dtc = , dtcvar = , fitc =, pitc = {
              chunkSize <-  ceiling(maxMemory/model$k)}
            )
    
    if (!("S" %in% names(model)))
      varsigma <-  matrix(0, dim(X)[1], model$d) 
    else
      varsigma <-  matrix(0, dim(X)[1], 1) 
    
    startVal <-  1 
    endVal <-  chunkSize 
    if (endVal> dim(X)[1])
      endVal <- dim(X)[1] 
    
    
    while (startVal <= dim(X)[1])
    {
      indices <- startVal:endVal 
      
# % Compute kernel for new point.
      switch (EXPR = model$approx, 
              ftc = {
                KX_star <-  kernCompute(model$kern, model$X, X[indices, ])}, 
              dtc = , dtcvar = , fitc = , pitc = {
                KX_star <-  kernCompute(model$kern, model$X_u, X[indices, ])}   
              )
      
      
# % Compute variances.
      if ((!("isSpherical" %in% names(model))) || model$isSpherical)
      {# % Compute diagonal of kernel for new point.
        diagK <-  kernDiagCompute(model$kern, X[indices, ]) 
        switch (EXPR = model$approx,
                ftc = {
                  Kinvk <-  model$invK_uu%*%KX_star }, 
                dtc =, dtcvar =, fitc =, pitc = {
                  Kinvk <-  (model$invK_uu - (1/model$beta)*model$Ainv)%*%KX_star}
                )
        varsig <-  diagK - (colSums(KX_star*Kinvk)) 
        if ("beta" %in% names(model))
          varsig <-  varsig + (1/model$beta) 
        
        if (!("S" %in% names(model)))
          varsigma[indices, ] <-  repmat(varsig, 1, model$d) 
        else
          varsigma[indices, ] <-  varsig 
        
      } else {
        diagK <-  kernDiagCompute(model$kern, X[indices, ]) 
        for (i in 1:model$d)
        {
          cat("gpPosteriorVar.R \n")
          ind <-  model$indexPresent[[i]] 
          switch (EXPR = model$approx,
                  ftc = {
                    Kinvk <- model$invK_uu[i]*KX_star[ind,]
                  },
                  stop("Non-spherical not yet implemented for any approximation other than 'ftc'") 
                  )
          varsigma[indices, i] <-  diagK - t(rowSums(KX_star[ind, ]*Kinvk)) 
        }
      }
      
      if (!("S" %in% names(model)))
      {
        varsigma[indices, ] <-  varsigma[indices, ]*repmat(model$scale*model$scale, length(indices), 1) 
      }
      
# % Prepare for the next chunk.
      startVal <-  endVal + 1 
      endVal <-  endVal + chunkSize 
      if (endVal > dim(X)[1])
        endVal <-  dim(X)[1] 
      
    }
    
    return (varsigma)
  }
