lvmLoadData <-
function(dataset, dataDir = system.file("extdata",package="vargplvm"))
{
# % LVMLOADDATA description not available.  
# % FORMAT
# % DESC 
# % description not available.

  #   baseDir <-system.file("extdata",package="vargplvm")
  switch (EXPR = dataset, 
          oil = {
            dataSetName <- paste(dataDir, "/3Class.mat", sep = "")
            data <-readMat(dataSetName)},
          mESC = {dataSetName <- paste(dataDir, "/mESC.mat", sep = "")
                  data <-readMat(dataSetName)})
  return(data)
}
