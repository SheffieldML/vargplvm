rbfard2KernGradX <-
function(kern, X, X2)
{
# % RBFARD2KERNGRADX Gradient of RBFARD2 kernel with respect to input locations.
# % FORMAT
# % DESC computes the gradident of the automatic relevance determination radial basis function
# % kernel with respect to the input positions where both the row
# % positions and column positions are provided separately.
# % ARG kern : kernel structure for which gradients are being
# % computed.
# % ARG x1 : row locations against which gradients are being computed.
# % ARG x2 : column locations against which gradients are being computed.
# % RETURN g : the returned gradients. The gradients are returned in
# % a matrix which is numData2 x numInputs x numData1. Where numData1 is
# % the number of data points in X1, numData2 is the number of data
# % points in X2 and numInputs is the number of input
# % dimensions in X.
# %
# % SEEALSO rbfard2KernParamInit, kernGradX, rbfard2KernDiagGradX
# %
# % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
# %
# % COPYRIGHT : Michalis K. Titsias, 2009
# 
# % KERN


gX <- array(0, dim = c(dim(X2)[1], dim(X2)[2], dim(X)[1])) 
for (i in 1:dim(X)[1]) 
{
  gX[, , i] <- rbfard2KernGradXpoint(kern, t(X[i,]), X2) 
}
return (gX)
}
