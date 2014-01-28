cmpndKernParamInit <-
function (kern, matlabway = FALSE) {
# % CMPNDKERNPARAMINIT description not available.  
# % FORMAT
# % DESC 
# % description not available.
  
  kern$nParams <- 0
  kern$transforms <- list()

  if ( ! ("comp" %in% names(kern)) )
    kern$comp <- list()

  for ( i in seq(along=kern$comp) ) {

    kern$comp[[i]] <- kernParamInit(kern$comp[[i]], matlabway = matlabway)
    kern$nParams <- kern$nParams + kern$comp[[i]]$nParams
    kern$comp[[i]]$index <- array()

    if ( "numBlocks" %in% names(kern$comp[[i]]) ) {
      if ( i==1 ) {
        kern$numBlocks <- kern$comp[[i]]$numBlocks
      } else {
        if ( (!("numBlocks" %in% names(kern))) | (kern$numBlocks!=kern$comp[[i]]$numBlocks) ) {
          stop("Compound of multi kernels with different numbers of blocks.")
        }
      }
    } else {
      if ( "numBlocks" %in% names(kern) )
        stop("Attempt to combine multi-kernel with non multi-kernel.")
    }
  }

  kern$paramGroups <- diag(1, nrow=kern$nParams, ncol=kern$nParams)

  kern$whiteVariance <- 0
  kern$isStationary <- TRUE

  for ( i in seq(along=kern$comp) ) {
    if ( !kern$comp[[i]]$isStationary )
      kern$isStationary <- FALSE

    if ( kern$comp[[i]]$type == "white" ) {
      kern$whiteVariance <- kern$whiteVariance + kern$comp[[i]]$variance
    } else {
      if ( "whiteVariance" %in% names(kern$comp[[i]]) ) {
        kern$whiteVariance <- kern$whiteVariance + kern$comp[[i]]$whiteVariance
      }
    }
  }

  return (kern)
  
}
