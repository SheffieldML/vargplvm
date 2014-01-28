linard2VardistPsi0Compute <- function (linard2kern, vardist)
{
# % LINARD2VARDISTPSI0COMPUTE  description not available.  
# % FORMAT
# % DESC 
# % description not available
  
  A <- linard2kern$inputScales 
  Psi0 <- sum(A*colSums((vardist$means*vardist$means) + vardist$covars)) 
  return (Psi0)
# %Psi0covs <- sum(A.*sum(vardist.covars,1))  
  # 
# %psi00 <- kernCompute(linard2kern,vardist.means) 
# %ok <- 0 
# %for n<-1:size(vardist.means,1)
# %    ok <- ok + trace(diag(A)*(vardist.means(n,:)'*vardist.means(n,:))) + trace(diag(A)*diag(vardist.covars(n,:)))  
# %end
  #                               %ok
  #                               %pause
  #                               %Psi0 <- Psi0means + Psi0covs 
}
