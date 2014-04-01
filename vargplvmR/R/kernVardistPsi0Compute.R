kernVardistPsi0Compute <-
function (kern, vardist)
{
# % KERNVARDISTPSI0COMPUTE description not available.  
# % FORMAT
# % DESC 
# % description not available

if (!(kern$type == "cmpnd"))
{
  fhandle <- paste(kern$type, "VardistPsi0Compute",sep="")
  Psi0 <- do.call(fhandle,list(kern, vardist))
}
else  #% the kernel is cmpnd
{
  fhandle <- paste(kern$comp[[1]]$type, "VardistPsi0Compute", sep ="")
  Psi0 <- do.call(fhandle, list(kern$comp[[1]], vardist))
  
  for (i in 2:length(kern$comp))
  {
      fhandle <- paste(kern$comp[[i]]$type, "VardistPsi0Compute", sep = "")
      Psi0 <- Psi0 + do.call(fhandle, list(kern$comp[[i]], vardist))
  }
}

return (Psi0)
}
