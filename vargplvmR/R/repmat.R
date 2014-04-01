repmat <-
function(X,m,n)
{
# % REPMAT R equivalent of repmat (matlab)
# % FORMAT
# % DESC 
# % description not available.

  if (is.matrix(X))
  {
    mx = dim(X)[1]
    nx = dim(X)[2]
    out <- matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
  }
  else if (is.vector(X)) {
    mx = 1
    nx = length(X)
    out <- matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
  }
  else if (length(X) == 1)
  {
    out <- matrix(X,m, n)
  }
  return (out)
}
