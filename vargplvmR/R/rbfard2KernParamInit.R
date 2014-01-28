rbfard2KernParamInit <-
function (kern, matlabway = TRUE)
{
# % RBFARD2KERNPARAMINIT RBFARD2 kernel parameter initialisation.
# % The automatic relevance determination version of the radial basis
# % function kernel (RBFARD2) is a very smooth non-linear kernel and is a
# % popular choice for generic use.
# %
# % k(x_i, x_j) <- sigma2 * exp(-1/2 *(x_i - x_j)'*A*(x_i - x_j))
# %
# % The parameters are sigma2, the process variance (kern.variance), the
# % positive). 
# %
# % SEEALSO : rbfKernParamInit
# %
# % FORMAT
# % DESC initialises the automatic relevance determination radial basis function
# %  kernel structure with some default parameters.
# % ARG kern : the kernel structure which requires initialisation.
# % RETURN kern : the kernel structure with the default parameters placed in.
# %
# % SEEALSO : kernCreate, kernParamInit
# %
# % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
# %
# % COPYRIGHT : Michalis K. Titsias, 2009
# 
# % KERN


# % This parameter is restricted positive.
kern$variance <- 1 
kern$inputScales <- 0.999*matrix(1, 1, kern$inputDimension) 
kern$nParams <- 1 + kern$inputDimension 
kern$transforms[[1]] <- list()
kern$transforms[[1]]$index <- c(1:kern$nParams) 
kern$transforms[[1]]$type <- optimiDefaultConstraint("positive", matlabway = matlabway) 

kern$isStationary <- TRUE
return (kern)
}
