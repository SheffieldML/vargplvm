kernkernVardistPsi2Compute <-
function (kern1, kern2, vardist, Z)
{
# % KERNKERNVARDISTPSI2COMPUTE description not available.  
# % FORMAT
# % DESC 
# % description not available

if (kern2$type == kern1$type)
    stop("Kern1 and Kern2 cannot be the same")
else if ((!(kern1$type == "white")) && (!(kern2$type == "white")))
{
    swap <- 0 
    if (kern2$type == "rbfard2") 
    {
       swap <- 1 
    } else if ((kern2$type == "linard2") && (!(kern1$type == "rbfard2") ))
    {
      swap <- 1    
    }
    if (swap == 0)
    {
       fhandle <- paste(kern1$type, kern2$type, "VardistPsi2Compute", sep = "") 
       
       temp <- do.call(fhandle, list(kern1, kern2, vardist, Z))
       Psi2 <- temp$Psi2
       Pnovar <- temp$Pnobias
       Psi1 <- temp$Psi1
#     } else {
#        fhandle <- paste(kern2$type, kern1$type, "VardistPsi2Compute", sep = "") 
#        [Psi2, Pnovar, Psi1] <- do.call(fhandle, list(kern2, kern1, vardist, Z)) 
    }
}
else 
{
# the white kernel gives zeros
   Psi2 <- matrix(0,dim(Z)[1],dim(Z)[1])  
   Pnovar <- matrix(0, dim(vardist$means)[1], dim(vardist$means)[2]) 
   Psi1 <- matrix(0,dim(vardist$means)[1],dim(vardist$means)[2]) 
}
return (list(Psi2 = Psi2, Pnovar = Pnovar, Psi1 = Psi1))

}
