linard2biasVardistPsi2Gradient <- function (linardKern, biasKern, vardist, Z, covGrad)
{
# % LINARD2BIASVARDISTPSI2GRADIENT description not available.  
# % FORMAT
# % DESC 
# % description not available
  
# variational means
  N <- dim(vardist$means)[1] 
# inducing variables 
  M <- dim(Z)[1]  
  Q <- dim(Z)[2]
  
  #[Psi2, Pnobias, Psi1]
  l2bvp2c<- linard2biasVardistPsi2Compute(linardKern, biasKern, vardist, Z) 
  Psi2 <- l2bvp2c$Psi2
  Pnobias <- l2bvp2c$Pnobias
  Psi1 <- l2bvp2c$Psi1
  
# % inverse variances
  A <- linardKern$inputScales 
  
# % gradient for the bias parameter  
  gKern2 <- sum(Pnobias*covGrad) 
  Bnm <- biasKern$variance*matrix(1, dim(Psi1)[1], dim(Psi1)[2])  
  BPsi1Covg <- (Bnm%*%covGrad) 
  
  gKern <- rep(0, vardist$latentDimension)
  gVarmeans <- matrix(0, vardist$numData, vardist$latentDimension)
  gInd <- matrix(0, dim(Z)[1], vardist$latentDimension)
  for (q in 1:vardist$latentDimension)
  {
    gKern[q] <- sum((vardist$means[,q]%*%t(Z[,q]))*BPsi1Covg)
    
    gVarmeans[,q] <- A[q]*rowSums((matrix(1, vardist$numData,1)%*%t(Z[,q]))*BPsi1Covg) 
    gInd[,q] <- A[q]*t(colSums((vardist$means[,q]%*%matrix(1, 1,dim(Z)[1]))*BPsi1Covg)) 
    
  }
  
  gKern1 <- 2*matrix(gKern, nrow = 1)  
# % gVarmeans is N x Q matrix (N:number of data, Q:latent dimension)
# % this will unfold this matrix column-wise 
  gVarmeans <- 2*matrix(gVarmeans, nrow = 1) 
  
# % gInd is M x Q matrix (M:number of inducing variables, Q:latent dimension)
# % this will unfold this matrix column-wise 
  gInd <- 2*matrix(gInd, nrow = 1) 
  
  gVarcovars <- matrix(0, 1,prod(dim(vardist$covars)))
  
  return (list(gKern1 = gKern1, gKern2 = gKern2, gVarmeans = gVarmeans, 
               gVarcovars = gVarcovars, gInd = gInd))
}
