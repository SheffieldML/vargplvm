kernVardistPsi1Gradient <-
function (kern, vardist, Z, covGrad)
{
# % KERNVARDISTPSI1GRADIENT  description not available.  
# % FORMAT
# % DESC 
# % description not available

    if (!(kern$type == "cmpnd"))
    {
       fhandle <- paste(kern$type, "VardistPsi1Gradient", sep = "")
       vp1g <- do.call(fhandle, list(kern, vardist, Z, covGrad))    
#        [gKern, gVarmeans, gVarcovars, gInd]

#        % Transformations
       gKern <- paramTransformPsi1(kern, gKern)  
    } else { #% the kernel is cmpnd
       fhandle <- paste(kern$comp[[1]]$type, "VardistPsi1Gradient", sep = "") 
       vp1g <- do.call(fhandle, list(kern$comp[[1]], vardist, Z, covGrad)) 
       gKern <- vp1g$gKern
       gVarmeans <-vp1g$gVarmeans
       gVarcovars <- vp1g$gVarcovars
       gInd<- vp1g$gInd
#        % Transformations
       gKern <- paramTransformPsi1(kern$comp[[1]], gKern)  
       for (i in 2:length(kern$comp))
       {
           fhandle <- paste(kern$comp[[i]]$type, "VardistPsi1Gradient", sep = "") 
           vp1gi <- do.call(fhandle, list(kern$comp[[i]], vardist, Z, covGrad)) 
#            [gKerni, gVarmeansi, gVarcovarsi, gIndi]

#            % Transformations
           gKerni <- paramTransformPsi1(kern$comp[[i]], vp1gi$gKerni)  

           gVarmeans <- gVarmeans + vp1gi$gVarmeans
           gVarcovars <- gVarcovars + vp1gi$gVarcovars 
           gInd <- gInd + vp1gi$gInd 
           gKern <- c(gKern, vp1gi$gKern)  
       } 
    }

# % variational variances are positive (This should rather go to
# % vargplvmLogLikeGradients)
# %gVarcovars <- (gVarcovars(:).*vardist$covars(:))'  

# %-----------------------------------------------------
# % This applies transformations 
# % This must be done similar to kernGradient at some point 
# % function gKern <- paramTransformPsi1(kern, gKern)
# % %
# % % 
# % if strcmp(kern$type,'rbfard2') | strcmp(kern$type,'rbfardjit')
# %    gKern(1) <- gKern(1)*kern$variance 
# %    gKern(2:end) <- gKern(2:end).*kern$inputScales  
# % elseif strcmp(kern$type,'linard2')
# %    gKern(1:end) <- gKern(1:end).*kern$inputScales 
# % elseif strcmp(kern$type,'bias')
# %    gKern <- gKern*kern$variance    
# % else
# %    % do nothing
# % end
return (list(gKern = gKern, gVarmeans = gVarmeans, gVarcovars = gVarcovars, gInd = gInd))
}
