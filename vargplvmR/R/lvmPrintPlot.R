lvmPrintPlot <-
function (model, lbls, capName, experimentNo, colour = 0, savedir = getwd())
{
# % LVMPRINTPLOT Print latent space for learnt model.
# % FORMAT 
# % DESC prints a latent space repsresentation for an LVM model.
# % ARG model : the model to use for plotting the latent space.
# % ARG lbls : any lables that are available for plotting.
# % ARG capName : the name of the saved plots.
# % ARG experimentNo : the experiment number to assign to the files.
# % 
# % SEEALSO : lvmScatterPlot
# % 
# % COPYRIGHT : Neil D. Lawrence, 2006
  # 
# % MLTOOLS
#   model = mm
#   colour = 0
  
  #   cat("Do print here lvmPrintPlot.....\n")

  if (colour)
    lvmScatterPlotColor(model, lbls) 
  else
    lvmScatterPlot(model, lbls) 
  
  modelType <- model$type 
  substring(modelType[1], 1, 1) <- toupper(substring(modelType[1], 1, 1))
  substring(capName[1], 1, 1) <- toupper(substring(capName[1], 1, 1)) 
  fileName <-  paste("dem", capName, modelType, experimentNo, sep = "")
# % directory <- ['../tex/diagrams'] 
  #   directory = ['/home/jie/H-drive/Software/vargplvm/matlab/tex'];
  directory <-savedir
  printPlot(fileName, directory, directory) 
  
  # figure
  # clf
  # ax <- axes('position', [0.05 0.05 0.9 0.9]) 
  # hold on
  x11()
  par(bg="white")
  if (dim(model$X)[2]==2)
  {
    if (!is.null(lbls))
    {
      if (is.list(lbls) && !(lbls[[1]] == "connect"))
      {
        lvmTwoDPlot(model$X, lbls, getSymbols(dim(lbls[[1]])[2]), newPlot = TRUE) 
      }else if  (!((lbls == "connect")[1])) {
        lvmTwoDPlot(model$X, lbls, getSymbols(dim(lbls)[2]), newPlot = TRUE) 
      }
    } else {
      lvmTwoDPlot(model$X, lbls) 
    }
  } else if (dim(model$X)[2]==3) {
    cat("to do lvmPrintPlot lvmThreeDplot\n")
#     set(gca, 'cameraPosition', [-0.0901445 0.0899564 4.51826], ...
#         'CameraPositionMode', 'auto', ...
#         'CameraTarget', [-0.0901445 0.0899564 0.0139032], ...
#         'CameraTargetMode', 'auto', ...
#         'CameraUpVector', [0 1 0], ...
#         'CameraUpVectorMode', 'auto', ...
#         'CameraViewAngle', [6.60861], ...
#         'CameraViewAngleMode', 'auto') 
#     
#     hold on
#     if ~isempty(lbls) && ~strcmp(lbls, 'connect')
#     lvmThreeDPlot(model.X, lbls, getSymbols(size(lbls, 2))) 
#     else
#       lvmThreeDPlot(model.X, lbls) 
#     end
  }
  lim <- list()
  for (i in 1:dim(model$X)[2])
  {       
  x1Min <- min(model$X[, i]) 
  x1Max <- max(model$X[, i]) 
  x1Span <- x1Max - x1Min 
  x1Min <- x1Min - 0.05*x1Span 
  x1Max <- x1Max + 0.05*x1Span 
  lim[[i]] <- c(x1Min, x1Max) 
  }
  
#   set(ax, 'xLim', lim{1}) 
#   set(ax, 'yLim', lim{2}) 
  if ((dim(model$X)[2]) > 2)
  {
    cat("to do model$X[2] >2 lvmPrintPlot\n")
#     set(ax, 'zLim', lim{3}) 
  }
#   set(ax, 'fontname', 'arial') 
#   set(ax, 'fontsize', 20) 
  printPlot(paste(fileName, "NoGray", sep = ""), directory, directory) 
# %print('-depsc', ['../tex/diagrams/' fileName 'NoGray'])
# %print('-deps', ['../tex/diagrams/' fileName 'NoGrayNoColour'])
  
}
