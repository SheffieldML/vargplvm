rbfard2VardistPsi2Compute <-
function (rbfardKern, vardist, Z)
{
  # TO DO wither turn on parallel code
# % RBFARD2VARDISTPSI2COMPUTE one line description
# % FORMAT
# % DESC description
# % RETURN K : description
# % RETURN outKern : description
# % RETURN sumKern : description
# % RETURN Kgvar : description
# % ARG rbfardKern : the kernel structure associated with the rbfard2 kernel.
# % ARG vardist : description
# % ARG Z : description
# %
# %
# % SEEALSO : others
# %
# %
# % COPYRIGHT : Michalis K. Titsias, 2009
# %
# % COPYRIGHT : Neil D. Lawrence, 2009
# %
# 
# % VARGPLVM

# try
#     pool_open = matlabpool('size')>0 
# catch e
#     pool_open = 0 
# end

# % The conditions for the parallel code to run, is the workers pool to be
# % open, the parallel flag to be active and the number of datapoints N to be
# % larger than a reasonable threshold (otherwise there is unecessary
# % thread-communication overhead).
  # to do parallel
# if pool_open && (isfield(vardist,'parallel') && vardist.parallel) && size(vardist.means,1) > 15
#     [K, outKern, sumKern, Kgvar] = rbfard2VardistPsi2ComputePar(rbfardKern, vardist, Z) 
# else
    out <- rbfard2VardistPsi2ComputeOrig(rbfardKern, vardist, Z) 

    #return ( list(K = K, outKern = outKern, sumKern = sumKern, Kgvar = Kgvar))

return (out)
}
