whiteVardistPsi2Gradient <-
function (whitekern, vardist, Z, covGrad)
{
# % WHITEVARDISTPSI2GRADIENT Compute gradient of white variational PSI2.
# % FORMAT
# % DESC description here.
# % RETURN gKern :
# % RETURN gVarmeans :
# % RETURN gVarcovars :
# % RETURN gInd :
# % ARG whitekern : the kernel structure associated with the white kernel.
# % ARG vardist :
# % ARG Z : 
# % ARG covGrad : 
# %
# % SEEALSO : 
# %
# % COPYRIGHT : Michalis K. Titsias, 2009
# %

# % VARGPLVM
  
  gKern <- 0 
  gVarmeans <- matrix(0, 1, prod(dim(vardist$means)))  
  gInd <- matrix(0, 1, prod(dim(Z)))  
  gVarcovars <- matrix(0, 1,prod(dim(vardist$covars)))  

return (list(gKern = gKern, gVarmeans = gVarmeans, gVarcovars = gVarcovars, gInd = gInd))
}
