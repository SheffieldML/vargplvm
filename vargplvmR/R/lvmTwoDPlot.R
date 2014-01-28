lvmTwoDPlot <-
  function (X, lbl = NULL, symbol, newPlot = FALSE)
  {
# % LVMTWODPLOT Helper function for plotting the labels in 2-D.
# % FORMAT
# % DESC helper function for plotting an embedding in 2-D with symbols.
# % ARG X : the data to plot.
# % ARG lbl : the labels of the data point.
# % ARG symbol : the symbols to use for the different labels.
# %
# % SEEALSO : lvmScatterPlot, lvmVisualise
# %
# % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2008
    # 
# % MLTOOLS
    
    if(is.character(lbl) && lbl == "connect")
    {
      connect <- TRUE 
      lbl <- NULL 
    } else {
      connect <- FALSE 
    }
    
    if (is.list(lbl))
    {
      lblStr <- lbl[[2]] 
      lbl <- lbl[[1]]
      labelsString <- TRUE 
    } else {
      labelsString <- FALSE 
    }
    
    if (nargs() < 3)
    {
      if (is.null(lbl))
        symbol <- getSymbols(1) 
      else
        symbol <- getSymbols(dim(lbl)[2]) 
    }
    
    returnVal <- NULL 
    textReturnVal <- NULL 
    # nextPlot <- get(axisHand, 'nextplot') 
    for (i in 1:dim(X)[1])
    {
      # if (i == 2)
      # set(axisHand, 'nextplot', 'add') 
      
      if (!is.null(lbl))
        labelNo <- which(lbl[i, ] != 0) 
      else
        labelNo <- 1 
      
# %try 
      # returnVal <- c(returnVal,  plot(X[i, 1], X[i, 2], cex = 3, 
      #                                 col = symbol[[labelNo]][1], 
      #                                 pch = symbol[[labelNo]][2]))
      if (i == 1 && newPlot) {
        plot(X[i, 1], X[i, 2], cex = 3, col = symbol[[labelNo]][1], 
             pch = symbol[[labelNo]][2], xlim = c(min(X[,1]), max(X[,1])), 
             ylim = c(min(X[,2]), max(X[,2])))
      } else {
        points(X[i, 1], X[i, 2], cex = 3, col = symbol[[labelNo]][1], 
               pch = symbol[[labelNo]][2])
      }
      # to do
      # if (labelsString)
      # textReturnVal <- [textReturnVal  text(X(i, 1), X(i, 2), ['   ' lblStr{i}])] 
      # end
      # if (connect)
      # {
      # if (i>1)
      # line([X(i-1, 1) X(i, 1)], [X(i-1, 2) X(i, 2)]) 
      # }
    }
    
  }
