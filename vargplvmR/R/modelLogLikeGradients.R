modelLogLikeGradients <-
function (model)
{
# % MODELLOGLIKEGRADIENTS Compute a model's gradients wrt log likelihood.
# % FORMAT
# % DESC is a wrapper function to compute the gradients of the log
# % likelihood of a given model.
# % ARG model : the model for which likelihoods are computed.
# % RETURN g : teh gradients of the likelihood with respect to the
# % parameters.
# %
# % SEEALSO : modelCreate
# %
# % COPYRIGHT : Neil D. Lawrence, 2006, 2005
# 
# % MLTOOLS

fhandle <- paste(model$type, "LogLikeGradients", sep = "")

g <- do.call(fhandle, list(model))

if ("paramGroups" %in% names(model))
  g <- g*model$paramGroups

dim(g) <- c(1, length(g))
return (g)

}
