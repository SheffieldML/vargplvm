kernVardistPsi1Compute <-
function (kern, vardist, Z)
{
# % KERNVARDISTPSI1COMPUTE description not available.  
# % FORMAT
# % DESC 
# % description not available

if (!(kern$type == "cmpnd"))
{
  fhandle <- paste(kern$type, "VardistPsi1Compute", sep = "")
  temp <- do.call(fhandle,list(kern, vardist, Z))
  Psi1 <- temp$Psi1
  P <- temp$P
} else {#% the kernel is cmpnd 
  fhandle <- paste(kern$comp[[1]]$type, "VardistPsi1Compute", sep = "")
  temp <- do.call(fhandle, list(kern$comp[[1]], vardist, Z))
  Psi1 <- temp$K
  P <- temp$Knovar
  for (i in 2:length(kern$comp))
  {
      fhandle <- paste(kern$comp[[i]]$type, "VardistPsi1Compute", sep = "")
      temp <- do.call(fhandle, list(kern$comp[[i]], vardist, Z))
      Ptmp <- temp$K
      P <- temp$P
      Psi1 <- Psi1 + Ptmp
  }
}

return (list(Psi1 = Psi1, P =P ))
}
