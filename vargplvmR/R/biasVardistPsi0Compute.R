biasVardistPsi0Compute <-
function (biaskern, vardist)
{
# % BIASVARDISTPSI0COMPUTE one line description
# % FORMAT
# % DESC description
# % RETURN Psi0 : description
# % ARG biasKern : the kernel structure associated with the white kernel.
# % ARG vardist : description
# %
# %
# % SEEALSO : others
# %
# %
# % COPYRIGHT : Michalis K. Titsias, 2009
# %
# 
# % VARGPLVM

Psi0 <- vardist$numData*biaskern$variance 

return (Psi0)
}
