vardistCreate <-
function (X, Q, type, constraint)
{
# % VARDISTCREATE creates the structure of the variational distirbution over the latent values in the GP-LVM 
# % FORMAT
# % DESC
# % The variational distribution is assumed to be factorized over the 
# % the latent variables. Each factor each a Gaussain N(x_n|mu_n,S_n) 
# % with a diagonal covariance matrix  
# %
# % The structure of the vardist is similar to the structure of a kernel

if (nargs() == 3)
    constraint <- optimiDefaultConstraint("positive", matlabway = TRUE) 
else
    constraint <- "identity" 


vardist <- list()
vardist$type <- "vardist" 
vardist$vartype <- type  

vardist$numData <- dim(X)[1]  
vardist$latentDimension <- Q  
vardist$nParams <- 2*vardist$numData*vardist$latentDimension 
# %  number of training points
N <- dim(X)[1]

transforms <- list()
transforms[[1]] <-list()
transforms[[1]]$index <- c((N*Q+1):vardist$nParams) 
transforms[[1]]$type <- constraint 
vardist$transforms <- transforms

# % initialize the parameters
vardist$means  <- X 
# set.seed(setseed)
vardist$covars <- 0.1*matrix(1,N,Q) + 0.001*matrix(rnorm(N*Q),N,Q) 
vardist$covars[which(vardist$covars<0.05)] <- 0.05 
# %vardist$covars <- (eps^2)*ones(size(vardist$covars)) 
# %vardist$means <- randn(N,Q) 
# %pmeans <- randn(Q,N) 
# %pcovars <- randn(Q,N) 
# %params <- [pmeans(:)' pcovars(:)']  
# %vardist <- vardistExpandParam(vardist, params) 
return (vardist)
}
