whiteVardistPsi2Compute <-
function (whitekern, vardist, Z)
{
# % WHITEVARDISTPSI2COMPUTE one line description
# % FORMAT
# % DESC description
# % RETURN Psi2 : description
# % RETURN P : description
# % ARG whiteKern : the kernel structure associated with the white kernel.
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

Psi2 <- matrix(0, dim(Z)[1], dim(Z)[1]) 
P <- NULL
return (list(Psi2 = Psi2, P = P))
}
