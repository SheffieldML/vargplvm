gaussianPriorParamInit <-
function (prior)
{
# % GAUSSIANPRIORPARAMINIT Gaussian prior model's parameter initialisation.
# % FORMAT
# % DESC initialises the parameters of the Gaussian prior with some
# % default parameters.
# % ARG prior : prior structure to be initialised.
# % RETURN prior : prior structure with initial values in place.
# % 
# % SEEALSO : priorCreate, gammaPriorParamInit
# %
# % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
# 
# % PRIOR

prior$precision <- 1 

prior$transforms$index <- 1 
prior$transforms$type <- optimiDefaultConstraint("positive", 
matlabway = TRUE) 
prior$nParams <- 1 
return (prior)
}
