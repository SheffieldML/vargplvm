kernGradient <-
  function (kern, x, matlabway = FALSE, ...) {
# % KERNGRADIENT description not available.  
# % FORMAT
# % DESC 
# % description not available.
    funcName <- paste(kern$type, "KernGradient", sep="")
    func <- get(funcName, mode="function")
    g <- func(kern, x, matlabway = matlabway, ...)
    
    if (!matlabway)
    {      
      factors <- .kernFactors(kern, "gradfact")
      for (i in seq(along=factors))
        g[factors[[i]]$index] <- g[factors[[i]]$index]*factors[[i]]$val
    } else{      
      factors <- kernFactors(kern, "gradfact")
      g[factors$index] <- g[factors$index]*factors$val
    }
    return (g)
  }
