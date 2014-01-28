rbfard2VardistPsi2Gradient <-
function (rbfardKern, vardist, Z, covGrad)
{
# % RBFARD2VARDISTPSI2GRADIENT  description not available.  
# % FORMAT
# % DESC 
# % description not available
  
  # try
  #     pool_open <- matlabpool('size')>0 
  # catch e
  #     pool_open <- 0 
  # end
  
# % The conditions for the parallel code to run, is the workers pool to be
# % open, the parallel flag to be active and the number of datapoints N to be
# % larger than a reasonable threshold (otherwise there is unecessary
# % thread-communication overhead).
  # if pool_open && (isfield(vardist,'parallel') && vardist.parallel) && size(vardist.means,1) > 15
  #     [gKern, gVarmeans, gVarcovars, gInd] <- rbfard2VardistPsi2GradientPar(rbfardKern, vardist, Z, covGrad) 
  # else
  r2vp2g <- rbfard2VardistPsi2GradientOrig(rbfardKern, vardist, Z, covGrad)
  #     [gKern, gVarmeans, gVarcovars, gInd]
  # end
  
  return (r2vp2g)
}
