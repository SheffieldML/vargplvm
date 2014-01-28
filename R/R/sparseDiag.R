sparseDiag <-
function (d)
{
# % SPARSEDIAG Create a diagonal matrix that is sparse from a vector.
# % FORMAT
# % DESC creates a diagonal matrix that is sparse from a vector.
# % ARG d : the diagonal vector from which the sparse diagonal matrix
# % is formed.
# % RETURN D : the sparse diagonal matrix containing the vector as
# % its diagonal.
# %
# % SEEALSO : diag, spdiags
# %
# % COPYRIGHT : Neil D. Lawrence, 2005
  # 
# % NDLUTIL
  
  if (!is.vector(d))
  {
    if (is.matrix(d) && prod(dim(d)) != length(d))
      stop("Input must be a vector.")
  }
  D <- spdiags(d, 0, length(d), length(d))
# % % Can be made more efficient.
# % n = length(d);
# % D = spalloc(n, n, n);
# % for i = 1:n
# %   D(i, i) = d(i);
# % end
  
  return (D)
}
