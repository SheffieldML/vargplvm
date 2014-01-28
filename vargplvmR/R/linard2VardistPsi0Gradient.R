linard2VardistPsi0Gradient <- function (linard2Kern, vardist, covGrad)
{
# % LINARD2VARDISTPSI0GRADIENT description not available.  
# % FORMAT
# % DESC 
# % description not available
  
  A <-  linard2Kern$inputScales 
  gKern <-  covGrad*colSums((vardist$means*vardist$means) + vardist$covars)  
  
  gVarmeans <-  2*(vardist$means%*%(diag(A)))  
# %gVarmeans1 <-  2*(repmat(A,size(vardist.means,1),1).*vardist.means)  
  
  gVarcovars <-  matrix(1, dim(vardist$means)[1],1)%*%A  
  
  gVarmeans <-  covGrad*matrix(gVarmeans, nrow = 1)  
  gVarcovars <-  covGrad*matrix(gVarcovars, nrow = 1)
  return (list(gKern = gKern, gVarmeans = gVarmeans, gVarcovars = gVarcovars))
}
