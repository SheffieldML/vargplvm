rbfard2KernGradXpoint <-
function (kern, x, X2)
{
# % RBFARD2KERNGRADXPOINT Gradient with respect to one point of x.
# % FORMAT
# % DESC 
# % description not available


# print((diag(kern$inputScales)))
# cat("rbfard2KernGradXpoint")
# print (dim( (sqrt(diag(kern$inputScales)))))
scales <-(sqrt(diag(kern$inputScales))) 
gX <- matrix(0, dim(X2)[1], dim(X2)[2]) 
n2 <- dist2(X2%*%scales, x%*%scales) 
rbfPart <- kern$variance*exp(-n2*0.5) 
for (i in 1:dim(x)[2])
{
  gX[, i] <- as.vector(kern$inputScales[i]*(t(t(X2[, i])) - x[i])*rbfPart) 
}
return (gX)
}
