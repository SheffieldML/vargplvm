rbfard2VardistPsi0Gradient <-
function (rbfard2Kern, vardist, covGrad)
{
# % RBFARD2VARDISTPSI0GRADIENT  description not available.  
# % FORMAT
# % DESC 
# % description not available

  gKern <- matrix(0, 1, rbfard2Kern$nParams)  
  gKern[1] <- covGrad*vardist$numData 
  
  gVarmeans <- matrix(0, 1, prod(dim(vardist$means)))  
  gVarcovars <- matrix(0, 1, prod(dim(vardist$means)))  
  return(list(gKern = gKern, gVarmeans = gVarmeans, gVarcovars = gVarcovars))
}
