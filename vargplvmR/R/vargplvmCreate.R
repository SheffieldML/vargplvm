vargplvmCreate <-
  function (q, d, Y, options)
  {
    # # TO DO, ELSE # NOT CALLED
# % VARGPLVMCREATE Create a GPLVM model with inducing variables.
# % FORMAT
# % DESC creates a GP-LVM model with the possibility of using
# % inducing variables to speed up computation.
# % ARG q : dimensionality of latent space.
# % ARG d : dimensionality of data space.
# % ARG Y : the data to be modelled in design matrix format (as many
# % rows as there are data points).
# % ARG options : options structure as returned from
# % VARGPLVMOPTIONS. This structure determines the type of
# % approximations to be used (if any).
# % ARG enableDgtN: if set to true, the model will be created in D >> N
# % mode (if this is indeed the case).
# % RETURN model : the GP-LVM model.
# %
# % COPYRIGHT : Michalis K. Titsias and Neil D. Lawrence, 2009-2011
# % Modifications : Andreas C. Damianou, 2010-2011
# %
# % SEEALSO : vargplvmOptions
    # 
# % VARGPLVM
    
    # q = latentDim
    
    if (dim(Y)[2] != d)
      stop(paste("Input matrix Y does not have dimension ", d, sep=""))
    
    
    
    if ("enableDgtN"%in%names(options))
    { 
      enableDgtN <- options$enableDgtN 
    } else { 
      enableDgtN <- TRUE  #% default behaviour is to enable D >> N mode
    }
    
# % Datasets with dimensions larger than this number will run the
# % code which is more efficient for large D. This will happen only if N is
# % smaller than D, though.
    if (enableDgtN)
    {
      limitDimensions <- 5000  #% Default: 5000
    } else {
      limitDimensions <- 1e+10 
    }
    model <- list() 
    model$type <- "vargplvm" 
    model$approx <- options$approx 
    
    model$learnScales <- options$learnScales 
# %model.scaleTransform <- optimiDefaultConstraint('positive') 
    
    model$optimiseBeta <- options$optimiseBeta 
    
    if (("betaTransform"%in%names(options)) && is.list( options$betaTransform ))
    {    
      model$betaTransform <- options$betaTransform 
    } else { 
      model$betaTransform <-  optimiDefaultConstraint("positive", matlabway = TRUE) 
    }
    
    model$q <- q 
    model$d <- ncol(Y) 
    model$N <- nrow(Y) 
    
    model$optimiser <- options$optimiser 
    model$bias <- colMeans(Y) 
    model$scale <- rep(1, model$d) 
    
    if("scale2var1"%in%names(options))
    {
      if(options$scale2var1)
      {
        model$scale <- apply(Y, 2, sd) 
        model$scale[which(model$scale==0)] <- 1 
        if(model$learnScales)
          warning("Both learn scales and scale2var1 set for GP") 
        if("scaleVal"%in%names(options))
          warning("Both scale2var1 and scaleVal set for GP") 
      }
    }
    
    if("scaleVal"%in%names(options))
      model$scale <- repmat(options$scaleVal, 1, model$d) 
    
    
    model$y <- Y 
    model$m <- gpComputeM(model) 
    
    
    if (is.character(options$initX))
    {
      #     %%% The following should eventually be uncommented so as to initialize
      #     %%% in the dual space. This is much more efficient for large D.
      initFunc <- paste(options$initX, "Embed",sep = "") 
      X <- do.call(initFunc, list(model$m, q))$X
      
      # #NOT CALLED
    } else {
      if (dim(options$initX)[1] == dim(Y)[1] && dim(options$initX)[2] == q)
        X <- options$initX 
      else
        stop("options.initX not in recognisable form.") 
    }
    
    model$X <- X 
    
    
    model$learnBeta <- 1  #% Newly added: deafault value for learnBeta.
    
# % If the provided matrices are really big, the required quantities can be
# % computed externally (e.g. block by block) and be provided here (we still
# % need to add a switch to get that).
    # TO DO, NOT CALLED
    if (model$d > limitDimensions && model$N < limitDimensions)
    {
      cat("vargplvmCreate.R to do \n")
      #     model.DgtN = 1  % D greater than N mode on.
      #     fprintf(1, '# The dataset has a large number of dimensions (%d)! Switching to "large D" mode!\n',model.d) 
      # 
      #     % If we have test data, we can prepare some multiplications in
      #     % advance, obtain NxN matrices and then never store Y.
      # 
      
      #     % Keep the original m. It is needed for predictions.
      #     model.mOrig = model.m 
      # 
      #     % The following will eventually be uncommented.
      #     YYT = model.m * model.m'  % NxN
      #     % Replace data with the cholesky of Y*Y'.Same effect, since Y only appears as Y*Y'.
      #     %%% model.m = chol(YYT, 'lower')   %%% Put a switch here!!!!
      #     [U S V]=svd(YYT) 
      #     model.m=U*sqrt(abs(S)) 
      # 
      #     model.TrYY = sum(diag(YYT))  % scalar
      # 
    } else {
      #     % Trace(Y*Y) is a constant, so it can be calculated just once and stored
      #     % in memory to be used whenever needed.
      model$DgtN <- 0 
      model$TrYY <- sum(model$m * model$m) 
    }
    
    model$date <- Sys.Date() #%%%%
    
# %%% Also, if there is a test dataset, Yts, then when we need to take
# % model.m*my' (my=Yts_m(i,:)) in the pointLogLikeGradient, we can prepare in advance
# % A = model.m*Yts_m' (NxN) and select A(:,i).
# % Similarly, when I have m_new = [my m] and then m_new*m_new', I can find that
# % by taking m*m' (already computed as YY) and add one row on top: (m*my)'
# % and one column on the left: m*my.
# %
# % %%%%%%%%%%
    # 
# %%%% _NEW
    
    # TO DO
    if (is.list(options$kern))
      model$kern <- options$kern 
    else
      model$kern <- kernCreate(model$X, options$kern, matlabway = TRUE) 
    
# % check if parameters are to be optimised in model space (and constraints are handled
# % by optimiser)
    if (("notransform" %in% names(options)) && options$notransform == TRUE)
    {
      #     % store notransform option in model:
      #     % this is needed when computing gradients
      model$notransform <- true 
      model$vardist <- vardistCreate(X, q, "gaussian", "identity") 
    } else {
      model$vardist <- vardistCreate(X, q, "gaussian")  
    }
    
    if (options$approx == "dtcvar")
    {
      #         % Sub-sample inducing variables.
      model$k <- options$numActive 
      model$fixInducing <- options$fixInducing 
      if (options$fixInducing)
      {
        if (length(options$fixIndices)!=options$numActive)
          stop("Length of indices for fixed inducing variables must
                    match number of inducing variables")
            
        model$X_u <- model$X[options$fixIndices,] 
        model$inducingIndices <- options$fixIndices 
      } else {
        #             %%%NEW_: make it work even if k>N
        if (model$k <= model$N)
        {
          if (!("labels" %in% names(options)))
          {
#             set.seed(setseed)
            ind <- sample(1:model$N) 
            ind <- ind[1:model$k] 
            model$X_u <- model$X[ind,] 
          } else {
            cat("to do vargplvmCreate \n")
            # #                     % in the case that class labels are supplied, make sure that inducing inputs
            # #                     % from all classes are chosen
            #                     [idcs, nSmpls] = class_samples( options.labels, model.k ) 
            #                     
            #                     count = 1 
            #                     midx = [] 
            #                     for inds = idcs
            #                         ind   = inds{:} 
            #                         ind   = ind(randperm(numel(ind)))                             
            #                         idx  = ind(1:nSmpls(count)) 
            #                         
            #                         % test that there is no overlap between index sets
            #                         assert(isempty(intersect(midx, idx))) 
            #                         midx = [midx, idx] 
            #                         
            #                         count = count+1 
            #                     end    
            #                     model.X_u = model.X(midx,:) 
          }
        }  else {
          cat("to do vargplvmCreate 2 \n")
          #                 % TODO: sample from the variational distr. (this should probably go
          #                 % to the dynamics as well because the vardist. changes in the initialization for the dynamics.
          # 
          #                 %!!! The following code needs some more testing!
          #                 samplingInd=0  %% TEMP
          #                 if samplingInd
          #                     % This only works if k<= 2*N
          #                     model.X_u=zeros(model.k, model.q) 
          #                     ind = randperm(model.N) 
          #                     %ind = ind(1:model.N) 
          #                     model.X_u(1:model.N,:) = model.X(ind, :) 
          # 
          #                     % The remaining k-N points are sampled from the (k-N) first
          #                     % distributions of the variational distribution (this could be done
          #                     % randomly as well).
          #                     dif=model.k-model.N 
          #                     model.X_u(model.N+1:model.N+dif,:)= model.vardist.means(1:dif,:) + rand(size(model.vardist.means(1:dif,:))).*sqrt(model.vardist.covars(1:dif,:))   % Sampling from a Gaussian.
          #                 else
          #                     model.X_u=zeros(model.k, model.q) 
          #                     for i=1:model.k
          #                         %ind=randi([1 size(model.vardist.means,1)]) 
          #                         % Some versions do not have randi... do it with rendperm
          #                         % instead: 
          #                         % ceil(size(model.vardist.means,1).*rand) % alternative
          #                         ind=randperm(size(model.vardist.means,1)) 
          #                         ind=ind(1) 
          #                         model.X_u(i,:) = model.vardist.means(ind,:) 
          #                     end
          #                 end
        }
        #             %%%_NEW
      }
      model$beta <- options$beta 
    }
      
    if (is.list(options$prior)) 
    {
      model$prior <- options$prior 
    } else {
      if (!is.null(options$prior))
        model$prior <- priorCreate(options$prior) 
    }
    
    if (("notransform" %in% names(options)) && options$notransform == TRUE)
      model$prior$transforms$type <- "identity"     
    
# %model.vardist = vardistCreate(X, q, 'gaussian') 
    
    if (("tieParam" %in% names(options)) && !is.null(options$tieParam))
    {
      cat("to do vargplvmCreate 3\n")
      #   to do
      #     if (options$tieParam == "free")
      #         % paramsList =
      #     else
      #         startVal = model.vardist.latentDimension*model.vardist.numData + 1 
      #         endVal = model.vardist.latentDimension*model.vardist.numData 
      #         for q=1:model.vardist.latentDimension
      #             endVal = endVal + model.vardist.numData 
      #             index = startVal:endVal 
      #             paramsList{q} = index 
      #             startVal = endVal + 1 
      #         end
      #         model.vardist = modelTieParam(model.vardist, paramsList) 
      #     end
      #     %
    }
     
    return (model)
  }
