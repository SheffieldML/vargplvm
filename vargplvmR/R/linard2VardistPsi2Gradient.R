linard2VardistPsi2Gradient <- function (linard2Kern, vardist, Z, covGrad)
{
# % LINARD2VARDISTPSI2GRADIENT  description not available.  
# % FORMAT
# % DESC 
# % description not available
  
# % inverse variances
  A <- linard2Kern$inputScales
  
  Psi1 <- (linard2VardistPsi1Compute(linard2Kern, vardist, Z))$K 
  
  Amat <- diag(A) 
  K1 <- (t(vardist$means)%*%vardist$means + diag(colSums(vardist$covars))) 
  K <- Amat%*%K1%*%Amat 
  
# % TYPICAL WAY
# %sumS <- sum(vardist.covars,1) 
# %for q<-1:vardist.latentDimension
# %   sc <- sum(vardist.means(:,q)*ones(1,size(Z,1)).*Psi1,1)  
# %   
# %   ZZ <- Z(:,q)*Z(:,q)' 
  #       %      
  #       %   gKern(q)  <- sum(sum( (Z(:,q)*sc + (sumS(q)*A(q))*ZZ).*covGrad ))  
  #       %   
  #       %   gVarmeans(:,q) <-  A(q)*sum(Psi1*((ones(size(Z,1),1)*Z(:,q)').*covGrad),2) 
# %   
# %   gVarcovars(q) <- sum(sum((ZZ.*covGrad)))    
# %end
# %gKern <- 2*gKern(:)' 
# %gVarmeans <- 2*gVarmeans(:)'  
# %gVarcovars <- (A.^2).*gVarcovars 
# %gVarcovars <- ones(size(Psi1,1),1)*gVarcovars 
# %gVarcovars <- gVarcovars(:)' 
# %gInd <- covGrad*Z*K 
# %gInd <- 2*gInd(:)' 
# %%% end of typical way 
  
  #          % FAST WAY
  AZ <- Z%*%Amat 
  AZK <- AZ%*%K1 
  gKern <- colSums((covGrad%*%Z)*AZK) 
  gVarmeans <- (Psi1%*%(covGrad%*%Z))%*%Amat 
  AZ <- AZ%*%Amat 
  tmp <- colSums((covGrad%*%Z)*AZ)
  dim(tmp) <- c(1, length(tmp))
  gVarcovars <- matrix(1, dim(Psi1)[1],1)%*%tmp
  gInd <- covGrad%*%Z%*%K 
  
  gKern <- 2*matrix(gKern, nrow = 1) 
  gVarmeans <- 2*matrix(gVarmeans, nrow = 1)  
  gVarcovars <- matrix(gVarcovars, nrow = 1) 
  gInd <- 2*matrix(gInd, nrow = 1) 
# %%%  end of FAST WAY
  #          
  #          %sum(sum(abs(gKern1-gKern)))
  #          %sum(sum(abs(gVarmeans1 - gVarmeans)))
  #          %sum(sum(abs(gVarcovars1 - gVarcovars)))
  #          %sum(sum(abs(gInd1 - gInd)))
  #          %pause
  
  return (list(gKern = gKern, gVarmeans = gVarmeans, gVarcovars = gVarcovars, gInd = gInd))
}