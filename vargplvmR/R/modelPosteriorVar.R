modelPosteriorVar <-
function (model, X) 
{
# % MODELPOSTERIORVAR variances of the posterior at points given by X.
# % FORMAT
# % DESC returns the posterior  variance for a given set of
# % points.
# % ARG model : the model for which the posterior will be computed.
# % ARG x : the input positions for which the posterior will be
# % computed.
# % RETURN sigma : the variances of the posterior distributions.
# %
# % SEEALSO : modelCreate, modelPosteriorMeanVar
# %
# % COPYRIGHT : Neil D. Lawrence, 2009
  # 
# % MLTOOLS
  varExist <- FALSE 
  func <- paste(model$type, "PosteriorVar", sep = "") 
  if (exists(func, mode = "function"))
    varExist <- TRUE 
  
  if (!varExist)
  {
    cat("to do modelPosteriorVar \n")
    func <- paste(model$type, "PosteriorMeanVar", sep = "") 
  }
  fhandle <- func 
  
  # if str2num(version('-release'))>13
  # if varExist
  # varsigma <- fhandle(model, X) 
  # else
  #   [mu, varsigma] <- fhandle(model, X) 
  # end
  # else 
  if (varExist)
    varsigma <- do.call(fhandle, list(model, X)) 
  # else
  #   [mu, varsigma] <- feval(fhandle, model, X) 
  # end
  # end
  # end
  
  return (varsigma)
}
