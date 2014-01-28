printPlot <-
  function (fileName, directory, directoryHtml) 
  {
# % PRINTPLOT Print a plot to eps and png files.
# % FORMAT 
# % DESC prints a plot to the specified file name and directory.
# % ARG fileName : the base name of the file to use.
# % ARG directory : the directory to place the eps files in.
# % ARG directoryPng : the directory to place png the file in.
# % 
# % SEEALSO : preparePlot
# %
# % COPYRIGHT : Neil D. Lawrence, 2008
    # 
# % NDLUTIL
    
    # global printDiagram #bar <- function() {x <<- 1:5}
    
    if (is.null(printDiagram))
      printTrue <- TRUE # % User hasn't been explicit, assume print.
    else
      printTrue <- printDiagram  #% Use user preference.
    
    if (printTrue)
    {
      if (nargs() < 2) 
        directory <- "." 
      
      if (nargs() < 3)
        pngPlot <- FALSE 
      else
        pngPlot <- TRUE
      
      cat("Printing eps plot ...\n") 
      df <- dev.copy2eps(device = x11, file = paste(directory, "", fileName, ".eps", sep = "")) 
      
      if (pngPlot)
      {
        cat("Printing png plot ...\n") 
# % make smaller for PNG plot.
        #               dev.copy(device = x11)
#                 dev.copy(device = x11)#, file = paste(directoryHtml, fileName, 
        #                                             ".png", sep = ""), width = 480, height = 480)
        
        df <- dev.print(png, file=paste(directoryHtml, fileName, ".png", sep = ""), width = 700, height = 700)
      }
#       cat("Printing black and white eps plot ... not done yet\n") 
    }
  }
