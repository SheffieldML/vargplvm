cmpndKernDiagGradX <-
  function (kern, X, matlabway = FALSE) {
# % CMPNDKERNDIAGGRADX description not available.  
# % FORMAT
# % DESC 
# % description not available.

    X <- as.matrix(X)
    i <- 1
    
    if (!matlabway)
    {
      funcName <- paste(kern$comp[[i]]$type, "KernDiagGradX", sep="")
      func <- get(funcName, mode="function")
      
      if ( !is.na(kern$comp[[i]]$index) ) {
        gX <- array(0, dim=dim(X))
        gX[,kern$comp[[i]]$index,] <- func(kern$comp[[i]], X[,kern$comp[[i]]$index])
      } else {
        gX <- func(kern$comp[[i]], X)
      }
      
      for ( i in seq(2, length=(length(kern$comp)-1)) ) {
        if ( !is.na(kern$comp[[i]]$index) ) {
          gX[,kern$comp[[i]]$index] <- gX[,kern$comp[[i]]$index] +  func(kern$comp[[i]], X[,kern$comp[[i]]$index])
        } else {
          gX <- gX + func(kern$comp[[i]], X)
        }
      }
    } else {
      if ( !is.na(kern$comp[[i]]$index) ) {
        gX <- matrix(0, dim(X)[1], dim(X)[2])
        gX[,kern$comp[[i]]$index] <-
          kernDiagGradX(kern$comp[[i]], X[,kern$comp[[i]]$index])
      } else {
        gX <- kernDiagGradX(kern$comp[[i]], X)
      }
      
      for ( i in seq(2, length=(length(kern$comp)-1)) ) {
        if ( !is.na(kern$comp[[i]]$index) ) {
          gX[,kern$comp[[i]]$index] <- gX[,kern$comp[[i]]$index] + 
            kernDiagGradX(kern$comp[[i]], X[,kern$comp[[i]]$index])
        } else {
          gX <- gX + kernDiagGradX(kern$comp[[i]], X)
        }
      }
    }
    
    return (gX)
  }
