pca <-
function (data, N)
{
# % PCA Principal Components Analysis
# % FORMAT
# % DESC
# % PCCOEFF = PCA(DATA) computes the eigenvalues of the covariance
# % matrix of the dataset DATA and returns them as PCCOEFF.  These
# % coefficients give the variance of DATA along the corresponding
# % principal components.
# %
# % PCCOEFF = PCA(DATA, N) returns the largest N eigenvalues.
# %
# % [PCCOEFF, PCVEC] = PCA(DATA) returns the principal components as well
# % as the coefficients.  This is considerably more computationally
# % demanding than just computing the eigenvalues.
# %
# % SEE ALSO
# % EIGDEC, GTMINIT, PPCA
# %
# 
# % COPYRIGHT (c) Ian T Nabney (1996-2001)

if (missing(N)) 
   N <- dim(data)[2] 

if (N != round(N) || N < 1 || N > dim(data)[2])
   stop("Number of PCs must be integer, >0, < dim") 

# % Find the sorted eigenvalues of the data covariance matrix
PC <- eigdec(cov(data), N) 

# if evals_only
#    PCcoeff <- eigdec(cov(data), N) 
# else
#   [PCcoeff, PCvec] <- eigdec(cov(data), N) 

return (PC) 
}
