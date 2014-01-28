biasVardistPsi2Compute <-
function (biaskern, vardist, Z)
{
# % BIASVARDISTPSI2COMPUTE one line description
# % FORMAT
# % DESC description
# % RETURN Psi2 : description
# % RETURN P : description
# % ARG biasKern : the kernel structure associated with the bias kernel.
# % ARG vardist : description
# % ARG Z : description
# %
# %
# % SEEALSO : others
# %
# %
# % COPYRIGHT : Michalis K. Titsias, 2009
# %
# 
# % VARGPLVM


Psi2 <- repmat(vardist$numData*(biaskern$variance^2),dim(Z)[1],dim(Z)[1]) 
P <- NULL
return (list(Psi2 = Psi2, P = P))
}
