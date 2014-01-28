jitChol <-
  function (A, maxTries, only.values = TRUE)
  {
    #   to do
# % JITCHOL Do a Cholesky decomposition with jitter.
# % FORMAT
# % DESC attempts a Cholesky decomposition on the given matrix, if
# % matrix isn't positive definite the function gives a warning, adds
# % 'jitter' and tries again. At the first attempt the amount of
# % jitter added is 1e-6 times the mean of the diagonal. Thereafter
# % the amount of jitter is multiplied by 10 each time it is added
# % again. This is continued for a maximum of 10 times.
# % ARG A : the matrix for which the Cholesky decomposition is required.
# % ARG maxTries : the maximum number of times that jitter is added
# % before giving up (default 10).
# % RETURN U : the Cholesky decomposition for the matrix.
# %
# % FORMAT
# % DESC attempts a Cholesky decomposition on the given matrix, if
# % matrix isn't positive definite the function adds 'jitter' and tries
# % again. Thereafter the amount of jitter is multiplied by 10 each time
# % it is added again. This is continued for a maximum of 10 times.  The
# % amount of jitter added is returned.
# % ARG A : the matrix for which the Cholesky decomposition is required.
# % ARG maxTries : the maximum number of times that jitter is added
# % before giving up (default 10).
# % RETURN U : the Cholesky decomposition for the matrix.
# % RETURN jitter : the amount of jitter that was added to the
# % matrix.
# %
# % SEEALSO : chol, pdinv, logdet
# %
# % COPYRIGHT : Neil D. Lawrence, 2005, 2006
    # 
# % NDLUTIL
    
    if (nargs() < 2) #nargin < 2
      maxTries <- 10 
    nonPosDef <- 0
    jitter <- 0 
    for (i in 1:maxTries)
    {
      val = try({
        #     % Try --- need to check A is positive definite
        if (jitter == 0)
        {
          jitter <- abs(mean(diag(A)))*1e-6
          dig <- 10
          while(!isSymmetric(A) && dig > 4)
          {
            dig <- dig - 1
            A <- round(A, digits = dig)
          }
          #         val = try({
          UC <- chol(A)
          #           break}, silent = TRUE) 
          #         if (class(val) == "try-error")
          break
          #         else
          #           break
        } else {
          if (only.values)
            warning(paste("Matrix is not positive definite in jitChol, adding ", jitter, " jitter.", sep = ""))
          UC <- chol(Re(A+jitter*diag(dim(A)[1]))) 
          break
        }
      },silent = TRUE)
      if (class(val) == "try-error")
      {
#         UC <- chol(Re(A+jitter*diag(dim(A)[1]))) 
#         break
        # #     % Was the error due to not positive definite?
        #           nonPosDef <- 0 
        #           verString <- version 
        #           if (str2double(verString(1:3)) > 6.1)
        #           {
        #             [void, errid] <- lasterr 
        #             if (errid == "MATLAB:posdef")
                      nonPosDef <- 1 
        #           } else {
        #             errMsg <- lasterr 
        #             if findstr(errMsg, 'positive definite')
        #               nonPosDef <- 1 
        #           }
      }
      
      if (nonPosDef)
      {
        jitter <- jitter*10 
        if (i==maxTries)
          stop(paste("Matrix is non positive definite tried ",i,
                     " times adding jitter, but failed with jitter of ",
                     jitter, ". Increase max tries", sep =""))
        #       } else {
        #         stop(lasterr)
        }
    }
    if (only.values)
    {
      return (UC)
    } else {
      return (list(UC = UC, jitter = jitter))
    }
  }
