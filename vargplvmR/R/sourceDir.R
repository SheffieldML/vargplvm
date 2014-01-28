sourceDir <-
function(path = getwd(), trace = TRUE, ...) {
# % SOURCEDIR description not available.  
# % FORMAT
# % DESC 
# % description not available.
  for (nm in list.files(path, pattern = "\\.[Rr]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
