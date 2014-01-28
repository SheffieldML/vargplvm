cmpndKernExtractParam <-
  function (kern, only.values=TRUE,untransformed.values=FALSE,
            matlabway = TRUE) {
# % CMPNDKERNEXTRACTPARAM description not available.  
# % FORMAT
# % DESC 
# % description not available.

    startVal <- 1
    endVal <- 0
    
    if ( only.values ) {
      params <- c()
      
      for ( i in seq(along=kern$comp) ) 
        params <- c(params, kernExtractParam(kern$comp[[i]],
                                             untransformed.values=untransformed.values, 
                                             matlabway = matlabway))
      
    } else {
      storedTypes <- c()
      params <- c()
      paramNames <- c()
      origNames <- c()
      for ( i in seq(along=kern$comp) ) {
        paramsList <- kernExtractParam(kern$comp[[i]], only.values=only.values,
                                       untransformed.values=untransformed.values,
                                       matlabway = matlabway)
        params <- c(params, paramsList)
        kernName <- paste(kern$comp[[i]]$type, length(grep(kern$comp[[i]]$type, storedTypes))+1, sep="")
        paramName <- paste(kernName, names(paramsList), sep="_")
        origNames <- c(origNames, paramName)
        storedTypes <- c(storedTypes, kern$comp[[i]]$type)
      }
    }
    
    paramNames <- array()
    if ( "paramGroups" %in% names(kern) ) {
      paramGroups <- kern$paramGroups
      for ( i in seq(length.out=dim(paramGroups)[2]) ) {
        ind <- grep(1, paramGroups[,i])
        if ( !only.values ) {
          paramNames[i] <- origNames[ind[1]]
          for ( j in seq(2, length.out=length(ind)-1) )
            paramNames[i] <- paste(paramNames[i], origNames[ind[j]],sep="/")
        }
        
        paramGroups[ind[seq(2,length(ind),length=length(ind)-1)], i] <- 0
      }
    }
    
    params <- params%*%paramGroups
    if ( !only.values )
      names(params) <- paramNames
    
    return (params)
  }
