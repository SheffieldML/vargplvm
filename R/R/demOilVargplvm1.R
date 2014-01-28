vargplvm<-function(dataSetName ="oil",  
                   dataDir = system.file("extdata",package="vargplvm"), 
                   iters = 2, seed = 1e5, kernoptions = "rbfard2",
                   createDir = TRUE, runVargplvmDir = getwd(),
                   overwriteDir = FALSE)
{
# % DEMOILVARGPLVM1 Run variational GPLVM on oil data.
# % FORMAT
# % DESC 
# % description not available.

  if (createDir) {
    dirA<-newDir(runVargplvmDir = runVargplvmDir, overwriteFile = overwriteDir)
    dataDir <-dirA[2]
  } else {
    extdir<-system.file("extdata",package="batman")
    dirA<-c(extdir,extdir,extdir)
  }
  ctime <- format(Sys.time(), "%d_%b_%H_%M_%S")
  dirctime<-paste(dirA[3],"/",ctime,sep="")
  if(!file.exists(dirctime)) {
    dir.create(dirctime)
  }
  dir5<-paste(dirctime,"/",sep="")
  
  printDiagram <<- NULL
  setseed <- seed

  cwd <- getwd()
  setwd(dirctime)
#   iters <<-iters
  #   wr <-3
  #   cl<<-makeCluster(wr, type = "SOCK")
  #   registerDoSNOW(cl)
  
  # DEMOILVARGPLVM1 Run variational GPLVM on oil data.
  #tic
  # VARGPLVM
  
  # Fix seeds
  #randn('seed', 1e5) 
  #rand('seed', 1e5) 
  set.seed(setseed)
  # dataSetName <- "oil" 
  experimentNo <- 1 
  printDiagram <- 1 
  
  # load data
  
  data <- lvmLoadData(dataSetName, dataDir)
  Y <- data$DataTrn
  lbls <- data$DataTrnLbls
  
  # Set up model
  options <- vargplvmOptions("dtcvar") 
#   options$kern <- c("rbfard2", "bias", "white") # "linard2"
  options$kern <- c(kernoptions, "bias", "white")
  
  options$numActive <- 50  
  
  
  options$optimiser <- "scg" 
  latentDim <- 10 
  d <- ncol(Y) 
  
# % demo using the variational inference method for the gplvm model
  
  model <- vargplvmCreate(latentDim, d, Y, options) 
  
  model <- vargplvmParamInit(model, model$m, model$X) 

  model$vardist$covars <- 0.5*matrix(1, dim(model$vardist$covars)[1], dim(model$vardist$covars)[2]) +
    0.001*matrix(rnorm(dim(model$vardist$covars)[1]*dim(model$vardist$covars)[2]), 
                 dim(model$vardist$covars)[1], dim(model$vardist$covars)[2]) 
  model$learnBeta<-1 
  
# % Optimise the model.
  #iters <- 2 #%2000
  display <- 1 
#   stim<- system.time({
    model <- vargplvmOptimise(model, display, iters) 
#   })[3]
#   cat("vargplvmOptimise ")
#   print(stim)
  
  
  capName <- dataSetName 
  substring(capName[1], 1, 1) <- toupper(substring(capName[1], 1, 1)) 
  modelType <- model$type 
  substring(modelType[1], 1, 1) <- toupper(substring(modelType[1], 1, 1)) 
  save(model, file = paste(dir5, "dem", capName, modelType, experimentNo, ".RData", sep = "")) 
  
# % order wrt to the inputScales 
  
  mm <- vargplvmReduceModel(model,2) 
  
# %% plot the two largest twe latent dimensions 
  #   to do print
#exists("printDiagram", mode = "function") &&
  if ( printDiagram)
    lvmPrintPlot(mm, lbls, capName, experimentNo, savedir = dir5) 
  #   stopCluster(cl)
  if (file.exists("means.txt"))
  {
    file.remove("means.txt")
    file.remove("covars.txt")
    file.remove("asPlus1.txt")
    file.remove("aDasPlus1.txt")
    file.remove("ZmZm.txt")
    file.remove("covGrad.txt")
    
    file.remove("partInd2.txt")
    file.remove("partA2.txt")
    file.remove("gVarmeans.txt")
    file.remove("gVarcovars.txt")    
  }

  setwd(cwd)
  return (mm)
  # ts <- toc
}
