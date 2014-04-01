whiteKernParamInit <-
  function (kern, matlabway = FALSE)
  {
# % WHITEKERNPARAMINIT WHITE kernel parameter initialisation.
# % The white noise kernel arises from assuming independent Gaussian
# % noise for each point in the function. The variance of the noise is
# % given by the kern.variance parameter. The variance parameter is
# % constrained to be positive, either by an exponential
# % transformation (default) or, if the flag use_sigmoidab is set, 
# % by a sigmoid transformation with a customizable output range.
# % 
# % This kernel is not intended to be used independently, it is provided
# % so that it may be combined with other kernels in a compound kernel.
# %
# % SEEALSO : cmpndKernParamInit
# %
# % FORMAT
# % DESC initialises the white noise
# %  kernel structure with some default parameters.
# % ARG kern : the kernel structure which requires initialisation.
# % RETURN kern : the kernel structure with the default parameters placed in.
# %
# % SEEALSO : kernCreate, kernParamInit
# %
# % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
# %
# % COPYRIGHT : Jaakko Peltonen, 2011
# %
# % COPYRIGHT : Antti Honkela, 2012
    # 
# % KERN
    
    
    kern$variance <- exp(-2) 
    kern$nParams <- 1 
    
    if (!matlabway)
    {
      kern$paramNames <- c("variance")
      kern$transforms <- list(list(index=c(1), type="positive"))
    } else {
# % The white-noise kernel can be set to use a ranged sigmoid
# % (sigmoidab) transformation for the variance, instead of a plain
# % exponential transformation.
      kern$transforms[[1]]<- list()
      if (("options" %in% names(kern)) && 
        ("use_sigmoidab" %in% names(kern$options)) && 
        (kern$options$use_sigmoidab==1))
      {
        kern$transforms[[1]]$index <- 1 
        kern$transforms[[1]]$type <- "sigmoidab" # %optimiDefaultConstraint("positive") 
        kern$transforms[[1]]$transformsettings <- c(0, 1e6)
      } else if (("options" %in% names(kern)) && 
        ("boundedParam" %in% names(kern$options)) && kern$options$boundedParam)
      {
        kern$transforms[[1]]$index <- 1 
        kern$transforms[[1]]$type <- optimiDefaultConstraint("bounded", matlabway = matlabway) 
        kern$transforms[[1]]$transformsettings <- c(0, 1e6) 
      } else {
        kern$transforms[[1]]$index <- 1 
        kern$transforms[[1]]$type <- optimiDefaultConstraint("positive", matlabway = matlabway) 
      }  
    }
    kern$isStationary <- TRUE 
    return (kern)
  }
