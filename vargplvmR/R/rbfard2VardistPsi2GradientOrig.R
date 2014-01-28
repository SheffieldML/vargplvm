rbfard2VardistPsi2GradientOrig <-
  function (rbfardKern, vardist, Z, covGrad)
  {
# % RBFARD2VARDISTPSI2GRADIENTORIG description not available.  
# % FORMAT
# % DESC 
# % description not available

    N <- dim(vardist$means)[1] 
# %  inducing variables 
    M <- dim(Z)[1]  
    Q <- dim(Z)[2]
    
# % evaluate the kernel matrix 
    r2vp2c <- rbfard2VardistPsi2Compute(rbfardKern, vardist, Z) 
    # [K, outKern, sumKern, Kgvar]
    
# % inverse variances
    A <- rbfardKern$inputScales 
    
# % gradient wrt variance of the kernel 
    gKernvar <- 2*sum(r2vp2c$Kgvar*covGrad)   
    
    
# % 1) line compute 0.5*(z_mq + z_m'q) for any q and store the result in a "M x Q x M" 
# %  matrix where M is the number of inducing points and Q the latent dimension
# % 2) line compute the z_mq - z_m'q, for any q
    ZmZm  <- array(0, dim = c(M,Q,M)) 
    ZmDZm <- array(0, dim = c(M,Q,M)) 
    for (q in 1:dim(Z)[2])
    {
      ZmZm[,q,] <- 0.5*(array(Z[,q], dim=c(length(Z[,q]), 1, M)) + array(rep(Z[,q],rep(M, length(Z[,q]))),dim=c(M, 1, M))) 
      ZmDZm[,q,] <- array(Z[,q], dim = c(length(Z[,q]), 1, M)) - array(rep(Z[,q],rep(M, length(Z[,q]))),dim=c(M, 1, M))
    }
    
# % compute the terms 2 a_q s_nq^2 + 1, for n and q and srore the result in a 
# % "N x Q" matrix
    asPlus1 <- 2*(repmat(A, N, 1)*vardist$covars) + 1  
# % compute the terms a_q/(2 a_q s_nq^2 + 1), for n and q and store the result in a 
# % "N x Q" matrix
    aDasPlus1 <- repmat(A, N, 1)/asPlus1  
    
    covGrad <- (rbfardKern$variance^2)*(covGrad*r2vp2c$outKern) 
    covGrad <- array(covGrad, dim = c(M, 1, M)) 
    sumKern <- array(r2vp2c$sumKern, dim = c(M, 1, M)) 
    Amq <- repmat(A, M, 1) 
    
    prod1 <- sumKern*covGrad
    temp <- array(0, dim = dim(ZmDZm))
    for (q in 1:dim(temp)[3])
    {
      temp[,,q] <- matrix(prod1[,,q], length(prod1[,,q]), Q) 
    }
    partInd1 <- - Amq*apply((ZmDZm*temp), 1:2, sum) #sum(ZmDZm.*repmat(sumKern*covGrad,[1 Q 1]),3) 
    partInd2 <- matrix(0, M, Q) 
    
    partA1 <- - 0.25*colSums(apply((ZmDZm*ZmDZm*temp),1:2, sum)) 
    partA2 <- matrix(0, 1, Q) 
    
    gVarcovars <- matrix(0, N, Q)  
    gVarmeans <- matrix(0, N, Q) 
    
# % Compute the gradient wrt lengthscales, variational means and variational variances  
# % For loop over training points  
    
    write.table(vardist$means,file="means.txt", sep = "\t", row.names=FALSE,col.names=TRUE,quote=FALSE)
    write.table(vardist$covars,file="covars.txt", sep = "\t", row.names=FALSE,col.names=TRUE,quote=FALSE)
    write.table(asPlus1,file="asPlus1.txt", sep = "\t", row.names=FALSE,col.names=TRUE,quote=FALSE)
    write.table(aDasPlus1,file="aDasPlus1.txt", sep = "\t", row.names=FALSE,col.names=TRUE,quote=FALSE)
    write.table(ZmZm,file="ZmZm.txt", sep = "\t", row.names=FALSE,col.names=TRUE,quote=FALSE)
    write.table(covGrad,file="covGrad.txt", sep = "\t", row.names=FALSE,col.names=TRUE,quote=FALSE)

    out <-.Call("vargplvm", as.integer(M), as.integer(N), as.integer(Q), A, 
                # vardist$means, vardist$covars, asPlus1, aDasPlus1,
                as.integer(dim(ZmZm)[2]), as.integer(dim(covGrad)[2]), PACKAGE = "vargplvm")

    
    partInd2<-read.table("partInd2.txt", header=FALSE,sep="\t")
    partInd2 <- t(partInd2)
    partA2<-read.table("partA2.txt", header=FALSE,sep="\t")
    partA2 <- t(partA2)
    gVarmeans<-read.table("gVarmeans.txt", header=FALSE,sep="\t")
    gVarmeans <- t(gVarmeans)
    gVarcovars<-read.table("gVarcovars.txt", header=FALSE,sep="\t")
    gVarcovars <- t(gVarcovars)
        
    gInd <- partInd1 + 2*partInd2  
    
    gKernlengcs <- partA1 - partA2  
    gKern <- c(gKernvar, gKernlengcs) 
    
# % gVarmeans is N x Q matrix (N:number of data, Q:latent dimension)
# % this will unfold this matrix column-wise 
# %gVarmeans <- gVarmeans'  
    gVarmeans <- matrix(gVarmeans, nrow = 1)  
    
# % gVarcovars is N x Q matrix (N:number of data, Q:latent dimension)
# % this will unfold this matrix column-wise 
# %gVarcovars <- gVarcovars'  
    #   print(dim)
    gVarcovars <- matrix(gVarcovars, nrow = 1) 
    
# % gInd is M x Q matrix (M:number of inducing variables, Q:latent dimension)
# % this will unfold this matrix column-wise 
# %gInd <- gInd'  
    gInd <- matrix(gInd, nrow = 1)  
    return (list(gKern = gKern, gVarmeans = gVarmeans, gVarcovars = gVarcovars, gInd = gInd))
  }
