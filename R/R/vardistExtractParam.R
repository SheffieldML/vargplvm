vardistExtractParam <-
function (vardist, only.values=TRUE, untransformed.values = FALSE)
{
# % VARDISTEXTRACTPARAM Extract a parameter vector from a vardist structure.
# % FORMAT
# % DESC extracts a parameter vector from a given VARDIST structure.
# % ARG model : the model from which parameters are to be extracted.
# % RETURN params : the parameter vector extracted from the model.
# %
# % DESC does the same as above, but also returns parameter names.
# % ARG model2 : the model structure containing the information about
# % the model.
# % RETURN params : a vector of parameters from the model.
# % RETURN names : cell array of parameter names.
# %
# % COPYRIGHT : Neil D. Lawrence, 2005, 2006
# %
# % SEEALSO : vardistCreate, vardistExpandParam, modelExtractParam
# 
# % VARGPLVM



# %means = vardist.means'  
# %covs = vardist.covars' 

# % the variational means and diagonal covariances obtained COLUMN-WISE 
params <- c(matrix(vardist$means, 1, prod(dim(vardist$means))), 
           matrix(vardist$covars, 1, prod(dim(vardist$covars))))
# To do
# % names
if (!only.values) #nargout > 1  
  cat("name vardistExtractParam\n")
#     for i=1:size(vardist$means,1)
#     for j=1:size(vardist$means,2)
#         varmeanNames{i,j} <- ['varmean(' num2str(i) ', ' num2str(j) ')'] 
#     end
#     end
#     for i<-1:size(vardist$means,1)
#     for j<-1:size(vardist$means,2)
#         varcovNames{i,j} <- ['varcov(' num2str(i) ', ' num2str(j) ')'] 
#     end
#     end
#     names(params) <- c(matrix(varmeanNames, ncol = 1), matrix(varcovNames, ncol= 1)) 
# end

# % Check if parameters are being optimised in a transformed space.
if (!is.null(vardist$transforms))
{
  for (i in 1:length(vardist$transforms))
  {
    index <- vardist$transforms[[i]]$index 
    fhandle <- paste(vardist$transforms[[i]]$type, "Transform", sep = "")
#     fhandle <- vardist$transforms[[i]]$type$func
    params[index] <- do.call(fhandle, list(params[index], "xtoa", 
    matlabway = TRUE)) 
  }
}
return (params)
}
