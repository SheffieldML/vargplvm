whiteVardistPsi1Compute <-
function (whitekern, vardist, Z)
{
# % WHITEVARDISTPSI1COMPUTE one line description
# % FORMAT
# % DESC description
# % RETURN K : description
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

K <- matrix(0,dim(vardist$means)[1], dim(Z)[1])
P <- NULL
return (list(K = K, P = P))
}
