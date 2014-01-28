biasKernParamInit <-
  function (kern, matlabway = TRUE)
  {
# % BIASKERNPARAMINIT BIAS kernel parameter initialisation.
# % The bias kernel arises from placing a prior over the bias with a
# % variance given by the kern.variance parameter. It allows the
# % output function to move up and down. 
# %
# % This kernel is not intended to be used independently, it is
# % provided so that it may be combined with other kernels in a
# % compound kernel.
# %
# % SEEALSO : cmpndKernParamInit
# %
# % FORMAT
# % DESC initialises the bias
# %  kernel structure with some default parameters.
# % ARG kern : the kernel structure which requires initialisation.
# % RETURN kern : the kernel structure with the default parameters placed in.
# %
# % SEEALSO : kernCreate, kernParamInit
# %
# % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
# %
# % COPYRIGHT : Antti Honkela, 2012
    # 
# % KERN
    
    
    kern$variance <- exp(-2) 
    kern$nParams <- 1 
    kern$transforms[[1]] <- list()
    kern$transforms[[1]]$index <- 1 
    if (("options" %in% names(kern)) && 
      ("boundedParam" %in% names(kern$options)) && kern$options$boundedParam)
    {
      kern$transforms[[1]]$type <- optimiDefaultConstraint("bounded", matlabway = matlabway) 
      kern$transforms[[1]]$transformsettings <- c(0, 1e6) 
    } else {
      kern$transforms[[1]]$type <- optimiDefaultConstraint("positive", matlabway = matlabway) 
    }
    
    kern$isStationary <- TRUE
    return (kern)
  }
