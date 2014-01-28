cmpndKernExpandParam <-
  function (kern, params, matlabway = FALSE) {
# % CMPNDKERNEXPANDPARAM description not available.  
# % FORMAT
# % DESC 
# % description not available.

    if ( is.list(params) )
      params <- params$values
    params <- params %*% (t(kern$paramGroups))  ## Params still log-transformed at this point.
    startVal <- 1
    endVal <- 0
    kern$whiteVariance <- 0
    for ( i in seq(along=kern$comp) ) {
      endVal <- endVal+kern$comp[[i]]$nParams
      kern$comp[[i]] <- kernExpandParam(kern$comp[[i]], params[startVal:endVal], 
      matlabway = TRUE)
      startVal <- endVal+1
      if (!matlabway)
      {
        if ( "white" %in% kern$comp[[i]]$type ) {
          kern$whiteVariance <- kern$whiteVariance+kern$comp[[i]]$variance
        } else if ( "whiteVariance" %in% names(kern$comp[[i]]) ) {
          kern$whiteVariance <- kern$whiteVariance+kern$comp[[i]]$whiteVariance
        }
      } else {
        if ("white" == substr(kern$comp[[i]]$type, 1, 5)) {
          kern$whiteVariance <- kern$whiteVariance+kern$comp[[i]]$variance
        } else if ( "whiteVariance" %in% names(kern$comp[[i]]) ) {
          kern$whiteVariance <- kern$whiteVariance+kern$comp[[i]]$whiteVariance
        }
      }      
    }
    
    return (kern)
  }
