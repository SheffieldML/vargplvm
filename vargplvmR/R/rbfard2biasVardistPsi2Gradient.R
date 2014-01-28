rbfard2biasVardistPsi2Gradient <-
function (rbfardKern, biasKern, vardist, Z, covGrad)
{

# % RBFARD2BIASVARDISTPSI2GRADIENT  description not available.  
# % FORMAT
# % DESC 
# % description not available
  
# % variational means
N <- dim(vardist$means)[1] 
# %  inducing variables 
M <- dim(Z)[1]  
Q <- dim(Z)[2]  


# [Psi2, Pnobias, Psi1] 
r2bvp2c<- rbfard2biasVardistPsi2Compute(rbfardKern, biasKern, vardist, Z) 


# % inverse variances
A <- rbfardKern$inputScales 

# % gradient wrt variance of the rbfard2 kernel 
gKernvar <- sum(r2bvp2c$Psi2*covGrad)/rbfardKern$variance   
# % gradient for the bias parameter  
gKern2 <- sum(r2bvp2c$Pnobias*covGrad)  
Bnm <- biasKern$variance*matrix(1, dim(r2bvp2c$Psi1)[1], dim(r2bvp2c$Psi1)[2])  
BPsi1Covg <- r2bvp2c$Psi1*(Bnm%*%covGrad)  

gVarmeans <- matrix(0, dim(BPsi1Covg)[1], vardist$latentDimension)
gInd <- matrix(0, dim(BPsi1Covg)[2], vardist$latentDimension)
gKernlengcs <- rep(0, vardist$latentDimension)
gVarcovars <- gVarmeans

# % compute the gradient wrt lengthscales, variational means and variational variances  
for (q in 1:vardist$latentDimension)
{
    S_q <- t(t(vardist$covars[,q]))   
    Mu_q <- t(t(vardist$means[,q]))  
    Z_q <- t(Z[,q])  
    
#     % B3_q term (without A(q)  see report)
    B_q <- (repmat(Mu_q, 1, M) - repmat(Z_q, N, 1))/repmat((A[q]*S_q + 1), 1, M) 
    
#     % derivatives wrt variational means and inducing inputs 
    tmp <- B_q * BPsi1Covg
    
#     % variational means: you sum out the columns (see report)
    gVarmeans[,q] <- -A[q]*rowSums(tmp)  
    
#     % inducing inputs: you sum out the rows 
    gInd[,q] <- A[q]*t(colSums(tmp))
    
#     % 
    B_q <- (B_q*(repmat(Mu_q, 1,  M) - repmat(Z_q, N, 1))) 
    
#     % B1_q term (see report)
    B1_q <- (repmat(S_q, 1, M) + B_q)/repmat((A[q]*S_q + 1), 1, M) 
    
#     % gradients wrt kernel hyperparameters (lengthscales) 
    gKernlengcs[q] <- -0.5*sum(B1_q*BPsi1Covg)  
    
#     % gradient wrt variational covars (diagonal covariance matrices) 
    gVarcovars[,q] <- rowSums((BPsi1Covg/repmat((A[q]*S_q + 1), 1, M))*(A[q]*B_q - 1)) 
    
#     %
}
# %

gKern1 <- c(gKernvar, 2*gKernlengcs) 

# % gVarmeans is N x Q matrix (N:number of data, Q:latent dimension)
# % this will unfold this matrix column-wise
gVarmeans <- 2*matrix(gVarmeans, nrow = 1)  

# % gVarcovars is N x Q matrix (N:number of data, Q:latent dimension)
# % this will unfold this matrix column-wise 
gVarcovars <- repmat(A, N, 1)*gVarcovars 
gVarcovars <- matrix(gVarcovars, nrow = 1) 

# % gInd is M x Q matrix (M:number of inducing variables, Q:latent dimension)
# % this will unfold this matrix column-wise 
gInd <- 2*matrix(gInd, nrow = 1)  

return (list(gKern1 = gKern1, gKern2 = gKern2, gVarmeans = gVarmeans, 
             gVarcovars = gVarcovars, gInd =gInd))
}
