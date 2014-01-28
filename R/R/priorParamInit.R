priorParamInit <-
function (prior)
{
# % PRIORPARAMINIT Prior model's parameter initialisation.
# % FORMAT
# % DESC wrapper function for initialising prior distributions parameters.
# % ARG prior : structure to initialise.
# % RETURN prior : initialised prior structure.
# %
# % SEEALSO : priorCreate
# %
# % COPYRIGHT : Neil D. Lawrence, 2003
#   
# % PRIOR

prior <- do.call(paste(prior$type, "PriorParamInit",sep = ""), list (prior))
return (prior)
}
