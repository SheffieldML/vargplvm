rbfard2VardistPsi1Compute <-
function (rbfardKern, vardist, Z)
{
# % RBFARD2VARDISTPSI1COMPUTE  description not available.  
# % FORMAT
# % DESC 
# % description not available

# % variational means
N  <- dim(vardist$means)[1]
# %  inducing variables 
M <- dim(Z)[1] 

A <- rbfardKern$inputScales
         
argExp <- matrix(0, N, M) 
normfactor <- matrix(1, N, 1)
for (q in 1:vardist$latentDimension)
{
    S_q <- vardist$covars[,q]  
    normfactor <- normfactor*(A[q]*S_q + 1)
    Mu_q <- vardist$means[,q] 
    Z_q <- t(Z[,q])
    distan <- (repmat(matrix(Mu_q),1, M) - repmat(Z_q,N, 1))^2
    argExp <- argExp + repmat(matrix(A[q]/(A[q]*S_q + 1)), 1, M)*distan
}

normfactor <- normfactor^0.5
Knovar <- repmat(1/normfactor, 1, M)*exp(-0.5*argExp)
K <- rbfardKern$variance*Knovar

return (list(K = K, Knovar = Knovar, argExp = argExp))
}
