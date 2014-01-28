whiteKernExpandParam <-
  function (kern, params, matlabway = TRUE) {
# % WHITEKERNEXPANDPARAM description not available.  
# % FORMAT
# % DESC 
# % description not available.
    if ( is.list(params) )
      params <- params$values
    
    kern$variance <- params[1] ## linear domain param, i.e. the untransformed noise variance
    
    return (kern)
  }
