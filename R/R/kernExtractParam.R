kernExtractParam <-
  function (kern, only.values=TRUE, 
  untransformed.values=FALSE, matlabway = FALSE) {
# % KERNEXTRACTPARAM description not available.  
# % FORMAT
# % DESC 
# % description not available.
    funcName <- paste(kern$type, "KernExtractParam", sep="")
    func <- get(funcName, mode="function")
    
    params <- func(kern, only.values=only.values, 
    untransformed.values=untransformed.values,
    matlabway = matlabway)
    
    if ( any(is.nan(params)) )
      warning("Parameter has gone to NaN.")
    
    if (!matlabway)
    {
    if ( "transforms" %in% names(kern) && (length(kern$transforms) > 0)
         && !untransformed.values )
      for ( i in seq(along=kern$transforms) ) {
        index <- kern$transforms[[i]]$index
        funcName <- optimiDefaultConstraint(kern$transforms[[i]]$type)
        func <- get(funcName$func, mode="function")
        if (funcName$hasArgs)
          params[index] <- func(params[index], "xtoa", kern$transformArgs[[i]])
        else
          params[index] <- func(params[index], "xtoa")
      }
    } else {
    
    if ( "transforms" %in% names(kern) && (length(kern$transforms) > 0)
         && !untransformed.values )
      for ( i in seq(along=kern$transforms) ) {
        index <- kern$transforms[[i]]$index
        funcName <-  paste(kern$transforms[[i]]$type, "Transform", sep = "")
        func <- get(funcName, mode="function")
        if (("transformsettings" %in% names(kern$transforms[[i]])) && 
          (length(kern$transforms[[i]]$transformsettings)>0))
          params[index] <- func(params[index], "xtoa", kern$transform[[i]]$transformsettings)
        else
          params[index] <- func(params[index], "xtoa", matlabway = matlabway)
      }

    }    
    return (params)
  }
