linard2VardistPsi1Compute <- function (linard2kern, vardist, Z)
{
# % LINARD2VARDISTPSI1COMPUTE  description not available.  
# % FORMAT
# % DESC 
# % description not available
  
  K <- kernCompute(linard2kern,vardist$means,Z);
  P <- NULL
  return (list(K = K, P = P))
}