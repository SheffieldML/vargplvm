biasVardistPsi2Gradient <-
function (biaskern, vardist, Z, covGrad)
{
# % BIASVARDISTPSI2GRADIENT Compute gradient of bias variational PSI2.
# % FORMAT
# % DESC description here.
# % RETURN gKern :
# % RETURN gVarmeans :
# % RETURN gVarcovars :
# % RETURN gInd :
# % ARG biaskern : the kernel structure associated with the bias kernel.
# % ARG vardist :
# % ARG Z : 
# % ARG covGrad : 
# %
# % COPYRIGHT : Michalis K. Titsias, 2009
# %

# % VARGPLVM

gKern <- (2*vardist$numData*biaskern$variance)*sum(matrix(1, dim(Z)[1],dim(Z)[1])*covGrad)  

gVarmeans <- matrix(0, 1, prod(dim(vardist$means)))  

gInd <- matrix(0, 1, prod(dim(Z)))  

gVarcovars <- matrix(0, 1, prod(dim(vardist$covars)))  

return (list(gKern = gKern, gVarmeans = gVarmeans, gVarcovars = gVarcovars, gInd = gInd))
}
