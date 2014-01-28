ppcaEmbed <-
function (Y, dims)
{
  ############## to do else # NOT CALLED, TO BE CONVERTED
# % PPCAEMBED Embed data set with probabilistic PCA.
# % FORMAT
# % DESC returns latent positions for a given data set via probabilistic
# % PCA.
# % ARG Y : the data set which you want the latent positions for.
# % ARG dims : the dimensionality of the latent space.
# % RETURN X : the latent positions.
# % RETURN sigma2 : the variance not explained by the latent positions.
# % RETURN W : the matrix required to invert the transformation, Y=X*W'.
# %
# % COPYRIGHT : Neil D. Lawrence, 2006
# %
# % SEEALSO : lleEmbed, isomapEmbed
  # 
# % MLTOOLS
  
  if (!any(is.nan(Y)))
  {
    if (dim(Y)[1]<dim(Y)[2])
    {
      Ymean <- colMeans(Y) 
      Ycentre <- matrix(0,dim(Y)[1],dim(Y)[2]) 
      for (i in 1:dim(Y)[1])
      {
        Ycentre[i, ] <- Y[i, ] -Ymean 
      }
      if (dim(Ycentre)[2]>30000)
      {
        #       % Bug in MATLAB 7.0 means you have to do this.
        innerY <- zeros(size(Ycentre, 1)) 
        
        for (i in 1:dim(Ycentre)[1])
          innerY[i, ] <- Ycentre[i, ]%*%t(Ycentre)
        
      } else {
        innerY <- Ycentre%*%t(Ycentre)
        #     
      }
      PC <- eigdec(innerY, dims)  
      v <- PC$values 
      u <- PC$vectors 
      v[which(v<0)] <- 0 
      X <- u[, 1:dims]*sqrt(dim(Y)[1]) 
      v <- v/sqrt(dim(Y)[1]) 
      sigma2 <- (sum(diag(innerY)) - sum(v))/(dim(Y)[2]-dims) 
      W <- t(X)%*%Ycentre 
    } else  {
      PCout <- pca(Y) 
      v <- PCout$values 
      u <- PCout$vectors 
      v[which(v<0)] <- 0 
      Ymean <- colMeans(Y) 
      Ycentre <- matrix(0,dim(Y)[1],dim(Y)[2]) 
      for (i in 1:dim(Y)[2])
        Ycentre[, i] <- Y[, i] - Ymean[i] 
      
      X <- Ycentre%*%u[, 1:dims]%*%diag(1/sqrt(v[1:dims])) 
      sigma2 <- mean(v[(dims+1):length(v)]) 
      W <- t(X)%*%Ycentre 
    }
  } 
  else { 
    cat("to do ppcaEmbed\n");
  # # NOT CALLED, TO BE CONVERTED
  # # % Hacky implementation of Probabilistic PCA for when there is missing data.
  #   iters <- 100 
  # # % Initialise W randomly
  #   d <- dim(Y)[2] 
  #   q <- dims 
  #   N <- dim(Y)[1] 
  #   W <- randn(d, q)*1e-3 
  #   sigma2 <- 1 
  #   mu <- matrix(0, d, 1) 
  #   for (i in 1:d) 
  #   {
  #     obs <- which(!is.nan(Y[, i])) 
  #     if (length(obs)>0)
  #       mu[i] <- mean(Y[obs, i]) 
  #     else
  #       mu[i] <- 0 
  #   }
  #   numObs <- sum(sum(~isnan(Y))) 
  #   for i <- 1:iters
  #   {
  #     M <- t(W)*W + sigma2*eye(q) 
  #     invM <- inv(M) 
  #     exp_xxT <- zeros(q) 
  #     exp_x <- matrix(0, N, q) 
  #     for (n in 1:N)
  #     {
  #       obs <- find(~isnan(Y(n, :))) 
  #       exp_x[n,] <- t(invM*t(W[obs,])*(t(Y[n, obs]) - mu[obs])) 
  #     }
  #     exp_xxT <- N*sigma2*invM + t(exp_x)*exp_x 
  #     s <- matrix(0, d, q) 
  #     s2 <- 0 
  #     for (n in 1:N)
  #     {
  #       obs <- which(!is.nan(Y[n,])) 
  #       subY <- matrix(0, dim(Y[n,])[2],dim(Y[n,])[1]) 
  #       subY[obs] <- t(Y[n, obs]) - mu[obs] 
  #       s <- s + (subY)*exp_x[n,] 
  #       s2 <- s2 + sum(t(Y[n, obs]) - mu[obs]).^2) - 2*exp_x[n,]*t(W[obs,])*(t(Y[n, obs]) - mu[obs]) 
  #     }
  #     W <- s*inv(exp_xxT) 
  #     sigma2 <- 1/(numObs)*(s2 + trace(exp_xxT*t(W)*W)) 
  #   }
#   X <- exp_x
  }
 
  out <- list() 
  out$X <- X 
  out$sigma2 <- sigma2 
  out$W <- W 
  return (out)
}
