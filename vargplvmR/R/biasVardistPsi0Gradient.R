biasVardistPsi0Gradient <-
function(biaskern, vardist, covGrad)
{
# % BIASVARDISTPSI0GRADIENT one line description
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
# %
# % SEEALSO : others
# %
# %
# % COPYRIGHT : Michalis K. Titsias, 2009
# %

# % VARGPLVM

gKern <- covGrad*vardist$numData 
 
gVarmeans <- matrix(0, 1, prod(dim(vardist$means)))  
gVarcovars <- matrix(0, 1, prod(dim(vardist$means)))  

return (list(gKern = gKern, gVarmeans = gVarmeans, gVarcovars = gVarcovars))
}
