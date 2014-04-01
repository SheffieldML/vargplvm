rbfard2VardistPsi1Gradient <-
function (rbfard2Kern, vardist, Z, covGrad)
{
# % RBFARD2VARDISTPSI1GRADIENT  description not available.  
# % FORMAT
# % DESC 
# % description not available
  
# % variational means
N <- dim(vardist$means)[1] 
# %  inducing variables 
M <- dim(Z)[1]  


# % evaluate the kernel matrix 
r2vp1c <- rbfard2VardistPsi1Compute(rbfard2Kern, vardist, Z) 
K_fu <- r2vp1c$K
Knovar <- r2vp1c$Knovar

# % inverse variances
A <- rbfard2Kern$inputScales 

# % gradient wrt variance of the kernel 
gKernvar <- sum(Knovar*covGrad)   

KfuCovGrad <- K_fu*covGrad 
gVarmeans <- matrix(0, dim(KfuCovGrad)[1], vardist$latentDimension)
gInd <- matrix(0, dim(KfuCovGrad)[2], vardist$latentDimension)
gKernlengcs <- rep(0, vardist$latentDimension)
gVarcovars <- gVarmeans
# % compute the gradient wrt lengthscales, variational means and variational variances  
for (q in 1:vardist$latentDimension)
{
    S_q <- t(t(vardist$covars[,q]))   
    Mu_q <- t(t(vardist$means[,q]))  
    Z_q <- t(Z[,q])  
    
#     % B3_q term (without A(q)  see report)
#     %B_q <- repmat(1./(A(q)*S_q + 1), [1 M]).*(repmat(Mu_q,[1 M]) - repmat(Z_q,[N 1])) 
    B_q <- (repmat(Mu_q, 1, M) - repmat(Z_q, N, 1))/repmat((A[q]*S_q + 1), 1, M) 
    
#     % derivatives wrt variational means and inducing inputs 
#     %tmp <- A(q)*((K_fu.*B_q).*covGrad) 
    tmp <- (B_q*KfuCovGrad) 
    
#     % variational means: you sum out the columns (see report)
    gVarmeans[,q] <- -A[q]*rowSums(tmp)  
    
#     % inducing inputs: you sum out the rows 
    gInd[,q] <- A[q]*t(colSums(tmp))  
    
#     % 
#     %B_q <- repmat(1./(A(q)*S_q + 1), [1 M]).*dist2(Mu_q, repmat(Z_q) 
    B_q <- (B_q*(repmat(Mu_q, 1, M) - repmat(Z_q, N, 1))) 
    
#     % B1_q term (see report)
#     %B1_q <- -(0.5./repmat((A(q)*S_q + 1), [1 M])).*(repmat(S_q, [1 M]) + B_q) 
    B1_q <- (repmat(S_q, 1, M) + B_q)/repmat((A[q]*S_q + 1), 1, M) 
    
#     % gradients wrt kernel hyperparameters (lengthscales) 
#     %gKernlengcs(q) <- sum(sum((K_fu.*B1_q).*covGrad))  
    gKernlengcs[q] <- -0.5*sum(B1_q*KfuCovGrad)  
    
#     % B2_q term (see report)  
#     %B1_q <- ((0.5*A(q))./repmat((A(q)*S_q + 1), [1 M])).*(B_q - 1)  
#   
#     %B1_q <- ((0.5*A(q))./repmat((A(q)*S_q + 1), [1 M])).*(A(q)*B_q - 1)  
#     
#     % gradient wrt variational covars (diagonal covariance matrices) 
#     %gVarcovars(:,q) <- sum((K_fu.*B1_q).*covGrad,2) 
    gVarcovars[,q] <- rowSums((KfuCovGrad/repmat((A[q]*S_q + 1), 1, M))*(A[q]*B_q - 1)) 
 
}

gKern <- c(gKernvar, gKernlengcs) 

# % gVarmeans is N x Q matrix (N:number of data, Q:latent dimension)
# % this will unfold this matrix column-wise 
# %gVarmeans <- gVarmeans'  
gVarmeans <- matrix(gVarmeans, nrow = 1)  

# % gVarcovars is N x Q matrix (N:number of data, Q:latent dimension)
# % this will unfold this matrix column-wise 
# %gVarcovars <- gVarcovars'  
gVarcovars <- 0.5*repmat(A, N, 1)*gVarcovars 
gVarcovars <- matrix(gVarcovars, nrow = 1) 

# % gInd is M x Q matrix (M:number of inducing variables, Q:latent dimension)
# % this will unfold this matrix column-wise 
# %gInd <- gInd'  
gInd <- matrix(gInd, nrow = 1 )
                                               
return (list(gKern = gKern, gVarmeans = gVarmeans, gVarcovars = gVarcovars, gInd = gInd))
}
