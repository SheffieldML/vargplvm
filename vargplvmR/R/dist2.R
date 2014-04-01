dist2 <-
function (x, c)
{
# % DIST2 Calculates squared distance between two sets of points.
# % FORMAT
# % DESC
# % D = DIST2(X, C) takes two matrices of vectors and calculates the
# % squared Euclidean distance between them.  Both matrices must be of
# % the same column dimension.  If X has M rows and N columns, and C has
# % L rows and N columns, then the result has M rows and L columns.  The
# % I, Jth entry is the  squared distance from the Ith row of X to the
# % Jth row of C.
# %
# % SEE ALSO
# % GMMACTIV, KMEANS, RBFFWD
# %
# 
# % COPYRIGHT (c) Ian T Nabney (1996-2001)
#   cat("x is")
#   print(dim(x))
#   print(class(x))
# x <- matrix(x, dim(x)[1], dim(x)[2])
# c <- matrix(c, dim(c)[1], dim(c)[2])
  ndata <- dim(x)
  ncentres <- dim(c)
if (ndata[2] != ncentres[2])
  stop("Data dimension does not match dimension of centres")


n2 <- t(matrix(1, ncentres[1], 1) %*% colSums(t(x^2))) + 
  matrix(1,ndata[1], 1) %*% colSums(t(c^2)) - 2*(x%*%t(c))

# % Rounding errors occasionally cause negative entries in n2
if (any(n2<0))
  n2[which(n2<0)] <- 0

return (n2)
}
