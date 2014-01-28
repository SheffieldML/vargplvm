kernFactors <-
function(kern, factorType)
{
# % KERNFACTORS Extract factors associated with transformed
# % optimisation space. "factorType" is one of 'atox', 'xtoa',
# % or 'gradfact', as follows.
# % FORMAT
# % DESC
# % 'atox': transform unrestricted input parameters to a suitably
# % restricted output range.
# %
# % 'xtoa': transform restricted parameters back to the unrestricted
# % space where e.g. gradient-based parameter optimization is done.
# %
# % 'gradfact': These factors are derivatives of the 
# % form dx/da, where "x" is the transformed parameter (for example
# % a parameter that is restricted to be positive by a suitable
# % transformation) and "a" is the untransformed parameter, which
# % usually can freely take any real values.
# %
# % COPYRIGHT : unknown original copyright, possibly Neil Lawrence
# %
# % COPYRIGHT : Jaakko Peltonen, 2011
# 
# % KERN
factors <- list()
factors$index <- NULL 
factors$val <- NULL 
if (length(kern$transforms) > 0)
{

  fhandle <- paste(kern$type, "KernExtractParam", sep ="") 
  params <- do.call(fhandle, list(kern)) 
  
# % Process each transformation used with this kernel. Each
# % transformation may affect several parameters.
  for (i in 1:length(kern$transforms))
  {
#     % Get the parameter indices involved in the i:th transformation
    index <- kern$transforms[[i]]$index 
    factors$index <- c(factors$index, index) 
    
#     % Get the handle of the transformation function of the i:th transformation
    fhandle <- paste(kern$transforms[[i]]$type, "Transform", sep = "") 
    
#     % If the transformation has been provided with specific
#     % settings (such as a custom output range), use the settings, 
#     % otherwise transform without settings
    if (("transformsettings" %in% names(kern$transforms[[i]])) && !is.null(kern$transforms[[i]]$transformsettings))
      factors$val <- c(factors$val,
                       do.call(fhandle, list(params[index], factorType, kern$transforms[[i]]$transformsettings))) 
    else
      factors$val <- c(factors$val, do.call(fhandle, list(params[index], factorType, matlabway = TRUE))) 
  }
    
}
return (factors)
}
