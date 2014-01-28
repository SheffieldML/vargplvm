kernVardistPsi0Gradient <-
function(kern, vardist, covGrad)
{
# % KERNVARDISTPSI0GRADIENT  description not available.  
# % FORMAT
# % DESC 
# % description not available
  
    if (!(kern$type == "cmpnd"))
    {
      cat("To do in kernVardistPsi0Gradient")
#        % 
#        fhandle <- str2func([kern$type 'VardistPsi0Gradient']) 
#        [gKern, gVarmeans, gVarcovars] <- fhandle(kern, vardist, covGrad)     
# 
#        % Transformations
#        gKern <- paramTransformPsi0(kern, gKern)  
#        %
    } else {#% the kernel is cmpnd
       fhandle <- paste(kern$comp[[1]]$type, "VardistPsi0Gradient", sep = "") 
#        [gKern, gVarmeans, gVarcovars] 
       vp0g <- do.call(fhandle, list(kern$comp[[1]], vardist, covGrad)) 
       gKern <- vp0g$gKern
       gVarmeans <- vp0g$gVarmeans
       gVarcovars <- vp0g$gVarcovars
#        % Transformations
       gKern <- paramTransformPsi0(kern$comp[[1]], gKern)  
#        %
       for (i in 2:length(kern$comp))
       {   
           fhandle <- paste(kern$comp[[i]]$type, "VardistPsi0Gradient", sep = "") 
#            [gKerni, gVarmeansi, gVarcovarsi] 
           kvp0gi <- do.call(fhandle, list(kern$comp[[i]], vardist, covGrad)) 

#            % Transformations
           gKerni <- paramTransformPsi0(kern$comp[[i]], kvp0gi$gKern)  

           gVarmeans <- gVarmeans + kvp0gi$gVarmeans  
           gVarcovars <- gVarcovars + kvp0gi$gVarcovars 
           gKern <- c(gKern, gKerni)  
#            %
      }
    } 
    
# % variational variances are positive (This should rather go to
# % vargplvmLogLikeGradients)
# %gVarcovars <- (gVarcovars(:).*vardist$covars(:))' 
return (list(gKern = gKern, gVarmeans = gVarmeans, gVarcovars = gVarcovars))
}
