lvmScatterPlot <-
function (model, YLbls, ax, dims, defaultVals) 
{
# % LVMSCATTERPLOT 2-D scatter plot of the latent points.
# % FORMAT
# % DESC produces a visualisation of the latent space with the given model.
# % ARG model : the model for which the scatter plot is being produced.
# % RETURN ax : the axes handle where the scatter plot was placed.
# %
# % DESC produces a visualisation of the latent space for the given model, 
# % using the provided labels to distinguish the latent points.
# % ARG model : the model for which the scatter plot is being produced.
# % ARG lbls : labels for each data point so that they may be given different
# % symbols. Useful when each data point is associated with a different
# % class.
# % RETURN ax : the axes handle where the scatter plot was placed.
# % 
# % DESC produces a visualisation of the latent space for the given model, 
# % using the provided labels to distinguish the latent points.
# % ARG model : the model for which the scatter plot is being produced.
# % ARG lbls : labels for each data point so that they may be given different
# % symbols. Useful when each data point is associated with a different
# % class.
# % ARG ax : the axes where the plot is to be placed.
# % RETURN ax : the axes handle where the scatter plot was placed.
# %
# % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
# %
# % SEEALSO : lvmVisualise, lvmTwoDPlot, lvmScatterPlotColor
  # 
# % MLTOOLS

  if (nargs() < 5)
  {
    defaultVals <- matrix(0, 1, dim(model$X)[2]) 
    if (nargs() < 4)
    {
      dims <- c(1, 2) 
      if (nargs()<3)
      {
        ax <- NULL 
        if (nargs() < 2)
        {
          YLbls <- NULL
        }
      }
    }
  }
  if (is.null(YLbls))
  {
    symbol <- getSymbols(1) 
  } else {
    if (is.list(YLbls))
      symbol <- getSymbols(dim(YLbls[1])[2]) 
    else
      symbol <- getSymbols(dim(YLbls)[2]) 
  }
  x1Min <- min(model$X[, dims[1]]) 
  x1Max <- max(model$X[, dims[1]]) 
  x1Span <- x1Max - x1Min 
  x1Min <- x1Min - 0.05*x1Span 
  x1Max <- x1Max + 0.05*x1Span 
  x1 <- seq(x1Min, x1Max, length.out = 150) 
  
  x2Min <- min(model$X[, dims[2]]) 
  x2Max <- max(model$X[, dims[2]]) 
  x2Span <- x2Max - x2Min 
  x2Min <- x2Min - 0.05*x2Span 
  x2Max <- x2Max + 0.05*x2Span 
  x2 <- seq(x2Min, x2Max, length.out = 150) 
  
# %if size(model.X, 2)==2
  
  # try
  # [X1, X2] <- meshgrid(x1, x2) 
  X1 <- repmat(x1,length(x2),1)
  X2 <- matrix(x2,length(x2), length(x1))
  XTest <- repmat(defaultVals, prod(dim(X1)), 1) 
  XTest[, dims[1]] <- c(X1) 
  XTest[, dims[2]] <- c(X2) 
  varsigma <- modelPosteriorVar(model, XTest) 
  posteriorVarDefined <- TRUE 
  # catch 
  # [lastMsg, lastId] <- lasterr 
  # disp(lastId)
  # if isoctave || strcmp(lastId, 'MATLAB:UndefinedFunction')
  # posteriorVarDefined <- false 
  # else
  #   rethrow(lasterror) 
  # end
  # end
  if (posteriorVarDefined)
  {
    d <- model$d 
    if (dim(varsigma)[2] == 1)
      dataMaxProb <- -0.5*d*log(varsigma) 
    else
      dataMaxProb <- -0.5*rowSums(log(varsigma)) 
    x11()
    
    C <- matrix(dataMaxProb, dim(X1)[1], dim(X2)[2]) 
    
# % Rescale it
    C <- C - (min(C)) 
    if (max(C) != 0)
    {
      C <- C/(max(C)) 
      C <- round(C*63) 
      if (x1[1] > x1[2])
      {
        x1 <- x1[length(x1):1]
        C <- C[dim(C)[1]:1,]
      }
      if (x2[1] > x2[2])
      {
        x2 <- x2[length(x2):1]
        C <- C[,dim(C)[2]:1]
      }
#       x11()
      image(x1, x2, t(C), col = gray.colors(256), cex.axis = 2)
#       x11()
#       image(x2, x1, C, col = gray.colors(256), cex.axis = 2)
    } #gray((0:32)/32))
    # 
# %[c, h] <- contourf(X1, X2, log10(reshape(1./varsigma(:, 1), size(X1))), 128)  
# % shading flat
    #   colormap gray 
# %colorbar
  }
  
  lvmTwoDPlot(model$X[, dims], YLbls, symbol) 
  switch (EXPR = model$type,
          dnet = {
            cat("switch lvmScatterPlot\n")
            # plot(model$X_u[, dims[1]], model$X_u[, dims[2]], 'g.')
          }) 
}
