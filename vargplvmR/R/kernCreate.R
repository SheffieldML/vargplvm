kernCreate <-
  function (X, kernelType, matlabway = FALSE)
  {
# % KERNCREATE Initialise a kernel structure.
# % FORMAT
# % DESC creates a kernel matrix structure given an design matrix of
# % DESC input points and a kernel type.
# % ARG X : Input data values (from which kernel will later be computed).
# % ARG type : Type of kernel to be created, some standard types are
# % 'lin', 'rbf', 'white', 'bias' and 'rbfard'. If a cell of the form
# % {'cmpnd', 'rbf', 'lin', 'white'} is used a compound kernel
# % based on the sum of the individual kernels will be created. The
# % 'cmpnd' element at the start of the sequence is
# % optional. Furthermore, {'tensor', 'rbf', 'lin'} can be used to
# % give a tensor product kernel, whose elements are the formed from
# % the products of the two indvidual kernel's elements and
# % {'multi', 'rbf', ...} can be used to create a block structured
# % kernel for use with multiple outputs.
# % Finally the form {'parametric', struct('opt1', {val1}), 'rbf'} can
# % be used to pass options to other kernels.
# % RETURN kern : The kernel structure.
# %
# % FORMAT
# % DESC creates a kernel matrix structure given the dimensions of
# % the design matrix and the kernel type.
# % ARG dim : input dimension of the design matrix (i.e. number of features
# % in the design matrix).
# % RETURN kern : The kernel structure.
# %
# % SEEALSO : kernParamInit
# %
# % COPYRIGHT: Neil D. Lawrence, 2006
# %
# % MODIFICATIONS: Antti Honkela, 2009, Andreas Damianou, 2012
    # 
    # KERN
    #   cat("debugging now ...\n")
    #   X <- model$X
    #   kernelType <- options$kern
    
    if (!matlabway)
    {
      if ( is.list(X) ) {
        dim <- array()
        for ( i in 1:length(X) ) {
          dim[i] <- dim(as.matrix(X[[i]]))[2]
          if ( (dim[i] == 1) & (dim(as.matrix(X[[i]]))[1] == 1) )
            dim[i] <- X[[i]]
        }
      } else {
        dim <- dim(as.matrix(X))[2]
        if ( (dim == 1) & (dim(as.matrix(X))[1] == 1) )
          dim <- X
      }
      
      if ( is.list(kernType) && kernType$type == "parametric" ) {
        kernOptions <- kernType$options
        kernType <- kernType$realType
      }
      
      if ( is.list(kernType) && ("options" %in% names(kernType)) ) {
        kernOptions <- kernType$options
      }
      
      if ( is.list(kernType) && ("complete" %in% names(kernType)) ) {
        if ( kernType$complete == 1 ) {
          kern <- kernType
        }
        
      } else if ( is.list(kernType) ) {
        
        kern <- list(inputDimension=dim, type=kernType$type)
        
        if (!is.null(kernOptions))
          kern$options <- kernOptions
        
        start <- 1    
        
        if ( kern$type == "multi" ) {
          for ( i in start:length(kernType$comp) ) {
            if ( is.list(kernType$comp) ) {
              iType <- kernType$comp[[i]]
            } else {
              iType <- kernType$comp[i]
            }
            
            if ( is.list(X) ) {
              kern$comp[[i-start+1]] <- kernCreate(X[[i-start+1]], iType)
              kern$diagBlockDim[i-start+1] <- dim(as.array(X[[i-start+1]]))[1]
            } else {
              kern$comp[[i-start+1]] <- kernCreate(X, iType)
            }
            
            kern$comp[[i-start+1]]$index = array()
          }
          
        } else if ( kern$type %in% c("cmpnd", "tensor", "translate",
                                     "selproj") )  {
          for ( i in start:length(kernType$comp) ) {
            if ( is.list(kernType$comp) ) {
              iType <- kernType$comp[[i]]
            } else {
              iType <- kernType$comp[i]
            }
            
            if (kern$type == "selproj") {
              if ( (dim(as.matrix(X))[2] == 1) && (dim(as.matrix(X))[1] == 1) )
                x_proj <- X-1
              else
                x_proj <- X[,-1]
              
              kern$comp[[i-start+1]] <- kernCreate(x_proj, iType)
            } else {
              kern$comp[[i-start+1]] <- kernCreate(X, iType)
            }
            kern$comp[[i-start+1]]$index = array()
          }
          
        } else if ( kern$type == "exp" ) {
          ## need double check
          if ( start == length(kernType$comp) ) {
            kern$argument <- kernCreate(X, kernType$comp[start])
          } else {
            kern$argument <- kernCreate(X, kernType$comp[start:length(kernType$comp)])
          }
        }
        
        kern <- kernParamInit(kern)
        
      } else {
        kern <- list(type=kernType, inputDimension=dim)
        
        if (!is.null(kernOptions))
          kern$options <- kernOptions
        
        kern <- kernParamInit(kern)
      }
      
      kern$Kstore <- matrix()
      kern$diagK <- matrix()      
      
      #   if (!is.null(kernOptions) && 'priors' %in% names(kernOptions)) {
      #     kern$priors <- list()
      #     for (k in seq_along(kernOptions$prior))
      #       kern$priors[[k]] <- priorCreate(kernOptions$prior[[k]])
      #   }
    } else {
      if (is.list(X))
      {
        dim <- array()
        for (i in 1:length(X)) {
          dim[i] <- dim(as.matrix(X[[i]]))[2]
          #     if ((dim[i] == 1) & (dim(as.matrix(X[[i]]))[1] == 
          #       1)) 
          #       dim[i] <- X[[i]]
        }
      } else {
# % This is a bit of a hack to allow the creation of a kernel without
# % providing an input data matrix (which is sometimes useful). If the X
# % structure is a 1x1 it is assumed that it is the dimension of the
# % input data.
        dim <- dim(X)[2] 
        if (dim == 1 && dim(X)[1] == 1) 
          dim <- X 
      }
      
      kern <- list()
      
      if (is.vector(kernelType) && (kernelType[1] == "parametric"))
      {
        kern$options <- kernelType[2] 
        kernelType <- kernelType[3]
      }
      
      if (is.vector(kernelType) && length(kernelType) > 1)
      {
        kern$inputDimension <- dim 
        switch (EXPR = kernelType[1],
                multi = {
# % multi output block based kernel.
                  start <- 2 
                  kern$type <- "multi"}, 
                tensor = {
# % tensor kernel type
                  start <- 2 
                  kern$type <- "tensor"}, 
                cmpnd = {
# % compound kernel type
                  start <- 2 
                  kern$type <- "cmpnd" },
                translate = {
# % translate kernel type
                  start <- 2 
                  kern$type <- "translate" },
                velotrans = {
# % velocity translate kernel type
                  start <- 2 
                  kern$type <- "velotrans" },
                exp = {
# % exponentiated kernel type
                  start <- 2 
                  kern$type <- "exp" },
                # otherwise
# % compound kernel type
{
  start <- 1 
  kern$type <- "cmpnd"}
                )
        
        switch (EXPR = kern$type,
                multi = {
                  kern$comp <- list()
                  for (i in start:length(kernelType))
                  {
                    if (is.list(X))
                    {                    
                      kern$comp[[i-start+1]] <- list()
                      kern$comp[[i-start+1]] <- kernCreate(X[[i-start+1]], kernelType[i], matlabway = matlabway) 
                      kern$diagBlockDim[[i-start+1]] <- dim(as.array(X[[i-start+1]]))[1] 
                    } else {
                      kern$comp[[i-start+1]] <- list()
                      kern$comp[[i-start+1]] <- kernCreate(X, kernelType[i], matlabway = matlabway) 
                    }
                    kern$comp[[i-start+1]]$index <- NULL
                  }
                }, 
                tensor =, cmpnd =, translate =, velotrans = {
# %%%---
                  #               cat("here \n")
                  if (kernelType[1] ==  "invcmpnd")
                  {
                    kern$type <- "invcmpnd" 
                    start <- 2 
                    if (is.list(X))
                    {
                      
# % Every kernel has its own input space, length(dim) should be
# % equal to length(kernelType(start:end))
                      if (length(dim) != length(kernelType[2:end]))
                        stop("For the invcmpnd kernel a separate input domain must be given for every compound") 
                      
                      
                      inds[[1]] <- 1:dim[1] 
                      kern$comp <- list()
                      kern$comp[[1]] <- list()
                      kern$comp[[1]] <- kernCreate(X[[1]], kernelType[2], matlabway = matlabway) 
                      kern$comp[[1]]$index <- inds[[1]] 
                      for (i in 2:length(kernelType)-1)
                      {
                        lastInd <- inds[[i-1]][length(inds[[i-1]])] 
                        inds[[i]] <- lastInd+1:lastInd + dim[i]  
                        kern$comp[[i]] <- list()
                        kern$comp[[i]] <- kernCreate(X[[i]], kernelType[i+1], matlabway = matlabway) 
                        kern$comp[[i]]$index <- inds[[i]] 
                      }
                    } else {
                      kern$comp <- list()
                      for (i in start:length(kernelType))                    
                      {                     
                        
                        kern$comp[[i-start+1]] <- list()
                        kern$comp[[i-start+1]] <- kernCreate(X, kernelType[i], matlabway = matlabway) 
                        kern$comp[[i-start+1]]$index <- NULL 
                      }
                    }
                  } else {
# %%%---
                    #                 cat("hre2\n")
                    kern$comp <- list()
                    for (i in start:length(kernelType))
                    {
                      kern$comp[[i-start+1]] <- list()
                      kern$comp[[i-start+1]] <- kernCreate(X, kernelType[i], matlabway = matlabway) 
                      kern$comp[[i-start+1]]$index <- NULL 
                    }
                  }},
                exp = {
                  if (start == length(kernelType))
                    kern$argument <- kernCreate(X, kernelType[start], matlabway = matlabway) 
                  else
                    kern$argument <- kernCreate(X, kernelType[start:end], matlabway = matlabway) 
                })
        kern <- kernParamInit(kern, matlabway = matlabway) 
      } else if (is.list(kernelType)) {
# % If a structure is passed, use it as the kernel.
        kern <- kernelType 
      } else {
        kern$type <- kernelType 
        if (is.list(X))
          kern$inputDimension <- dim[i] 
        else
          kern$inputDimension <- dim 
        
        kern <- kernParamInit(kern, matlabway = matlabway) 
      }
      kern$Kstore <- NULL 
      kern$diagK <- NULL
    }
    return (kern)
  }
