kernParamInit <-
function (kern, matlabway = FALSE) {
# % KERNPARAMINIT description not available.  
# % FORMAT
# % DESC 
# % description not available.
  funcName <- paste(kern$type, "KernParamInit", sep="")
  kern$transforms <- list()

  func <- get(funcName, mode="function")
  kern <- func(kern, matlabway = matlabway)  

  return (kern)
}
