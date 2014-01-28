cmpndKernGradient <-
  function (kern, x, matlabway = TRUE, x2, covGrad) {
# % CMPNDKERNGRADIENT description not available.  
# % FORMAT
# % DESC 
# % description not available.
    
    if ( missing(covGrad) ) 
    {
      covGrad <- x2
      misCovGrad <- TRUE
    }
    g <- array(0, dim(kern$paramGroups)[1])
    startVal <- 1
    endVal <- 0
    
    for ( i in seq(along=kern$comp) ) {
      endVal <- endVal + kern$comp[[i]]$nParams
      if ( !is.na(kern$comp[[i]]$index) ) {
        if (misCovGrad) {
          g[startVal:endVal] <- kernGradient(kern$comp[[i]], x[,kern$comp[[i]]$index], covGrad)
        } else {
          g[startVal:endVal] <- kernGradient(kern$comp[[i]], x[,kern$comp[[i]]$index], x2[,kern$comp[[i]]$index], covGrad)
        }
      } else {
        if (misCovGrad) {
          g[startVal:endVal] <- kernGradient(kern$comp[[i]], x, matlabway = matlabway, covGrad)
        } else {
          g[startVal:endVal] <- kernGradient(kern$comp[[i]], x, matlabway = matlabway, x2, covGrad)
        }
      }
      startVal <- endVal + 1       
    }
    
    g <- g %*% kern$paramGroups    
    
    return (g)
  }
