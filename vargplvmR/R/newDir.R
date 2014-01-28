newDir<-function(runVargplvmDir = getwd(), overwriteFile = FALSE)
{
# % NEWDIR Create a new directory with installed data files
# % FORMAT
# % DESC 
# % description not available.

  dirRoot<-paste(runVargplvmDir,"/runVargplvm",sep="")
  if(!file.exists(dirRoot)) {
    dir.create(dirRoot)
  }
  dirIP<-paste(dirRoot,"/VargplvmInput",sep="")
  if(!file.exists(dirIP)){
    dir.create(dirIP)
  }
  dirOP<-paste(dirRoot,"/VargplvmOutput",sep="")
  if(!file.exists(dirOP)) {
    dir.create(dirOP)
  }
  extdir<-system.file("extdata",package="vargplvm")
  dirbo<-paste(extdir,"/3Class.mat",sep="")
  cp<-file.copy(dirbo,dirIP, overwrite=overwriteFile)
#   dird<-paste(extdir,"/NMRdata.txt",sep="")
#   cp<-file.copy(dird,dirIP, overwrite=overwriteFile)
#   dirL<-paste(extdir,"/metabolitesList.csv",sep="")
#   cp<-file.copy(dirL,dirIP, overwrite=overwriteFile)
#   dirR<-paste(extdir,"/multi_data.csv",sep="")
#   cp<-file.copy(dirR,dirIP,overwrite=overwriteFile)
#   dirRU<-paste(extdir,"/multi_data_user.csv",sep="")
#   cp<-file.copy(dirRU,dirIP,overwrite=overwriteFile)
  
  dirA<-c(dirRoot, dirIP, dirOP)
  return (dirA)
}