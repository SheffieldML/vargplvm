whiteKernExtractParam <-
  function (kern, only.values=TRUE,
            untransformed.values=TRUE,
            matlabway = TRUE) {
# % WHITEKERNEXTRACTPARAM description not available.  
# % FORMAT
# % DESC 
# % description not available.
    params <- c(kern$variance)
    
    if ( !only.values ) {
      names(params) <- c("variance")
    }
    
    return (params)
  }
