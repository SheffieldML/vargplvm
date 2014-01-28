vargplvmObjective <-
function (params, model)
{

# % VARGPLVMOBJECTIVE Wrapper function for variational GP-LVM objective.
# % FORMAT
# % DESC provides a wrapper function for the variational GP-LVM, it
# % takes the negative of the log likelihood, feeding the parameters
# % correctly to the model.
# % ARG params : the parameters of the variational GP-LVM model.
# % ARG model : the model structure in which the parameters are to be
# % placed.
# % RETURN f : the negative of the log likelihood of the model.
# % 
# % SEEALSO : vargplvmCreate, vargplvmLogLikelihood, vargplvmExpandParam
# %
# % COPYRIGHT : Michalis K. Titsias, 2009-2011
# % COPYRIGHT : Neil D. Lawrence, 2009-2011
# 
# % VARGPLVM


model <- modelExpandParam(model, params)
f <- -vargplvmLogLikelihood(model)
# % fprintf(1,'# F: %.13f .\n',f); %%% DEBUG
return (f)
}
