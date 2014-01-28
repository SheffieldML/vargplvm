eigdec <-
function (x, N)
{
# % EIGDEC  Sorted eigendecomposition
# % FORMAT
# % DESC
# % EVALS = EIGDEC(X, N computes the largest N eigenvalues of the
# % matrix X in descending order.  [EVALS, EVEC] = EIGDEC(X, N) also
# % computes the corresponding eigenvectors.
# %
# % SEE ALSO
# % PCA, PPCA
# %
#                   
# % COPYRIGHT (c) Ian T Nabney (1996-2001)
                  
  if (N != round(N) || N < 1 || N > dim(x)[2])
  stop("Number of PCs must be integer, >0, < dim. \n")

  temp <- eigen(x) 
  out <- list() 
  out$values <- temp$values[1:N] 
  out$vectors <- temp$vectors[, 1:N] 
                  
  return (out)
}
