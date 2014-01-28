kernkernVardistPsi2Gradient <-
function (kern1, kern2, vardist, Z, covGrad)
{
# % KERNKERNVARDISTPSI2GRADIENT description not available.  
# % FORMAT
# % DESC 
# % description not available
if (kern2$type == kern1$type)
{
    stop("Kern1 and Kern2 cannot be the same")
    return  
} else if (!(kern1$type == "white") &&  !(kern2$type == "white"))
{
    swap <- 0 
    if (kern2$type == "rbfard2") 
       swap <- 1 
    else if ((kern2$type == "linard2") && !(kern1$type == "rbfard2")) 
       swap <- 1    

    if (swap == 0)
    {
       fhandle <- paste(kern1$type, kern2$type, "VardistPsi2Gradient", sep = "") 
       vp2g <- do.call(fhandle, list(kern1, kern2, vardist, Z, covGrad)) 
#        [gKern1, gKern2, gVarmeans, gVarcovars, gInd]
    } else {
       cat("Not checked, not runned here")
       fhandle <- paste(kern2$type, kern1$type, "VardistPsi2Gradient", sep = "") 
       vp2g <- do.call(fhandle, list(kern2, kern1, vardist, Z, covGrad)) 
#        [gKern1, gKern2, gVarmeans, gVarcovars, gInd]
    }
#     % !!! Do not take transformation here: they are applied in the kernkernVardistPsi2Gradient
#     %gKern1 <- paramTransformPsi2(kern1, gKern1) 
#     %gKern2 <- paramTransformPsi2(kern2, gKern2) 
#     % variational variances are positive  
#     %gVarcovars <- (gVarcovars(:).*vardist$covars(:))' 
#    %
} else {
#     % the white kernel gives zeros
    gKern1 <- matrix(0, 1, kern1$nParams)  
    gKern2 <- matrix(0, 1, kern2$nParams)  
    gVarmeans <- matrix(0, 1, prod(dim(vardist$means)))  
    gInd <- matrix(0, 1, prod(dim(Z)))  
    gVarcovars <- matrix(0, 1, prod(dim(vardist$covars))) 
    vp2g <-list(gKern1 = gKern1, gKern2 = gKern2,
                gVarmeans = gVarmeans,
                gVarcovars = gVarcovars, gInd =gInd)
#    %
}
return (vp2g)
}
