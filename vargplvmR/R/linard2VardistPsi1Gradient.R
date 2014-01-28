linard2VardistPsi1Gradient <- function (linard2Kern, vardist, Z, covGrad)
{
# % LINARD2VARDISTPSI1GRADIENT  description not available.  
# % FORMAT
# % DESC 
# % description not available

#   linard2Kern = kern$comp[[1]]
  
  A <- linard2Kern$inputScales 
  dim(A) <- c(1, length(A))
# %
# % TYPICAL WAY
# %for q<-1:vardist.latentDimension
# %   % 
# %   gKern(q) <- sum(sum((vardist.means(:,q)*(Z(:,q)')).*covGrad)) 
# %   
# %   gVarmeans(:,q) <- A(q)*sum((ones(vardist.numData,1)*Z(:,q)').*covGrad,2) 
# %   gInd(:,q) <- A(q)*sum((vardist.means(:,q)*ones(1,size(Z,1))).*covGrad,1)' 
# %   %
# %end
# %%% end of typical way            
# % FAST WAY
  AA <- matrix(1, dim(vardist$means)[1], 1)%*%A  
  covVarm <- t(covGrad)%*%vardist$means 
  gKern <- colSums(covVarm*Z)  
  gVarmeans <- AA*(covGrad%*%Z) 
  AA <- matrix(1, dim(Z)[1],1)%*%A 
  gInd <- AA*covVarm 
  
# %sum(sum(abs(gKern1-gKern)))
# %sum(sum(abs(gVarmeans1 - gVarmeans)))
# %sum(sum(abs(gInd1 - gInd)))
# %pause
  
  gKern <- matrix(gKern, nrow = 1)   
# % gVarmeans is N x Q matrix (N:number of data, Q:latent dimension)
# % this will unfold this matrix column-wise 
  gVarmeans <- matrix(gVarmeans, nrow = 1)
  
# % gInd is M x Q matrix (M:number of inducing variables, Q:latent dimension)
# % this will unfold this matrix column-wise 
  gInd <- matrix(gInd, nrow = 1)  
  
  gVarcovars <- matrix(0, 1, prod(dim(vardist$covars)))  
  return (list(gKern = gKern, gVarmeans = gVarmeans, gVarcovars = gVarcovars, gInd = gInd))
}