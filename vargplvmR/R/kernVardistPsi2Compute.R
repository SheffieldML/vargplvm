kernVardistPsi2Compute <-
  function (kern, vardist, Z)
  {
    
# % KERNVARDISTPSI2COMPUTE description not available.  
# % FORMAT
# % DESC    
# % compute first the "square" terms of Psi2 
#     print(kern$type)
#     cat("\n")
    if (!(kern$type == "cmpnd"))
    {   
      fhandle <- paste(kern$type, "VardistPsi2Compute", sep = "")
      temp <- do.call(fhandle, list(kern, vardist, Z))
#         cat("HERE ")
#         print(fhandle)
      Psi2 <- temp$Psi2
      P <- temp$P
    } else {#% the kernel is cmpnd
      fhandle <- paste(kern$comp[[1]]$type, "VardistPsi2Compute", sep = "")
      temp <- do.call(fhandle, list(kern$comp[[1]], vardist, Z))
      #   cat("HERE ")
      #   print(fhandle)
      if (kern$comp[[1]]$type == "rbfard2")
      {
        Psi2 <- temp$K
        P <- temp$outKern
      } else {
        Psi2 <- temp$Psi2
        P <- temp$P
      }
      for (i in 2:length(kern$comp))
      {
        fhandle <- paste(kern$comp[[i]]$type, "VardistPsi2Compute", sep="")
        temp <- do.call(fhandle, list(kern$comp[[i]], vardist, Z))
        Ptmp <- temp$Psi2
        P <- temp$P
        Psi2 <- Psi2 + Ptmp 
      }
    }
    
# % compute the cross-kernel terms of Psi2
    if ((kern$type == "cmpnd") && (length(kern$comp)>1))
    {
      for (i in 1:length(kern$comp))
      {
        j <- i + 1
        while (j <= length(kern$comp))
        {
          #        %[i j]      
          #       cat(length(kern$comp),"i = ", i, "j = ", j, "  ")
          temp <- kernkernVardistPsi2Compute(kern$comp[[i]], kern$comp[[j]], vardist, Z) 
          Ptmp <- temp$Psi2
          P <- temp$Pnovar
          Psi2 <- Psi2 + Ptmp
          j <- j + 1
        }
      }
    }
    return (list(Psi2 = Psi2, P = P))
  }
