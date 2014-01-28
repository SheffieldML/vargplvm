linard2VardistPsi2Compute <- function (linard2kern, vardist, Z)
{
# % LINARD2VARDISTPSI2COMPUTE  description not available.  
# % FORMAT
# % DESC 
# % description not available

  # 
  # 
# %Psi1 = linard2KernCompute(linard2kern, vardist.means, Z);
# %sqrtAS = sparse(diag(linard2kern.inputScales.*sqrt(sum(vardist.covars,1))));
# %Zsc = Z*sqrtAS;
# %P = Psi1'*Psi1;
# %Psi2 = P + Zsc*Zsc';
  
  ZA <- Z%*%(diag(linard2kern$inputScales))
  P <- t(vardist$means)%*%vardist$means + diag(colSums(vardist$covars))
  
  Psi2 <- ZA%*%P%*%t(ZA);
  return (list(Psi2 = Psi2, P = P))
}


