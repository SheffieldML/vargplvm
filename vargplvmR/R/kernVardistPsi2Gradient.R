kernVardistPsi2Gradient <-
function (kern, vardist, Z, covGrad)
{
# % KERNVARDISTPSI2GRADIENT description not available.  
# % FORMAT
# % DESC 
# % description not available

    if (!(kern$type == "cmpnd"))
    {
      cat("not implemented yet")
#       to do
#        % 
#        fhandle = str2func([kern$type 'VardistPsi2Gradient']) 
#        [gKern, gVarmeans, gVarcovars, gInd] = fhandle(kern, vardist, Z, covGrad)     
# 
#        % Transformations
#        gKern = paramTransformPsi2(kern, gKern)  
#        %
    } else {#% the kernel is cmpnd
       fhandle <- paste(kern$comp[[1]]$type, "VardistPsi2Gradient", sep = "") 
       vp2g <- do.call(fhandle, list(kern$comp[[1]], vardist, Z, covGrad)) 

       gKern <- vp2g$gKern
       gVarmeans <- vp2g$gVarmeans
       gVarcovars <- vp2g$gVarcovars
       gInd <- vp2g$gInd
#        [gKern, gVarmeans, gVarcovars, gInd]
#        % Transformations
       gKern <- paramTransformPsi2(kern$comp[[1]], vp2g$gKern)  

       for (i in 2:length(kern$comp))
       {
     
           fhandle <- paste(kern$comp[[i]]$type, "VardistPsi2Gradient", sep = "") 
#            stim <- system.time({
           vp2gi <- do.call(fhandle, list(kern$comp[[i]], vardist, Z, covGrad)) 
#            })[3]
#            cat(fhandle, " in kernVardistPsi2Gradient in for kern$comp ")
#            print(stim)
#          [gKerni, gVarmeansi, gVarcovarsi, gIndi] <-
#            % Transformations
           gKerni <- paramTransformPsi2(kern$comp[[i]], vp2gi$gKern)  
           gVarmeans <- gVarmeans + vp2gi$gVarmeans  
           gVarcovars <- gVarcovars + vp2gi$gVarcovars 
           gInd <- gInd + vp2gi$gInd 
           gKern <- c(gKern, gKerni)  
       }
    }


#     %% compute the cross-kernel terms of Psi2 
#     %
    index <- list()
    if ((kern$type == "cmpnd") && (length(kern$comp)>1))
    {
      index[[1]]<-1:kern$comp[[1]]$nParams  
      for (i in 2:length(kern$comp))
           index[[i]] <- (tail(index[[i-1]],1)+1):(tail(index[[i-1]],1)+kern$comp[[i]]$nParams) 

      for (i in 1:(length(kern$comp)-1))
      {
        for (j in (i+1):length(kern$comp))
        {

#            [gKerni, gKernj, gVarm, gVarc, gI] 
           kkvp2gi<- kernkernVardistPsi2Gradient(kern$comp[[i]], kern$comp[[j]], vardist, Z, covGrad) 
#            [gKern1, gKern2, gVarmeans, gVarcovars, gInd]
#            % Transformations
           gKerni <- paramTransformPsi2(kern$comp[[i]], kkvp2gi$gKern1) 
           gKernj <- paramTransformPsi2(kern$comp[[j]], kkvp2gi$gKern2) 
           gKern[index[[i]]] <- gKern[index[[i]]] + gKerni 
           gKern[index[[j]]] <- gKern[index[[j]]] + gKernj   

           gVarmeans <- gVarmeans + kkvp2gi$gVarmeans  
           gVarcovars <- gVarcovars + kkvp2gi$gVarcovars 
           gInd <- gInd + kkvp2gi$gInd 
#            %
        }
      }
    }

#     % variational variances are positive (This should rather go to
#                                           % vargplvmLogLikeGradients) 
# % gVarcovars <- (gVarcovars(:).*vardist$covars(:))' 
return (list(gKern = gKern, gVarmeans = gVarmeans, gVarcovars = gVarcovars, gInd = gInd))
}
