expTransform <-
  function (x, transform, matlabway = FALSE)
  {
# % EXPTRANSFORM Constrains a parameter to be positive through exponentiation.
# % FORMAT
# % DESC contains commands to constrain parameters to be positive via
# % exponentiation.
# % ARG x : input argument.
# % ARG y : return argument.
# % ARG transform : type of transform, 'atox' maps a value into
# % the transformed space (i.e. makes it positive). 'xtoa' maps the
# % parameter back from transformed space to the original
# % space. 'gradfact' gives the factor needed to correct gradients
# % with respect to the transformed parameter.
# % 
# % SEEALSO : negLogLogitTransform, sigmoidTransform
# %
# % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2007
    # 
# % OPTIMI

    if (!matlabway)
    {
      eps <- 2.2204e-16
      maxVal <- 4.3112e+15
      
      thre <- 36  ## threshold
      y <- array(0, dim(as.array(x)))
      
      if ( "atox" == transform ) {
        for ( ind in seq_along(as.array(x)) ) {
          if ( x[ind] > thre ) y[ind] <- maxVal else
            if ( x[ind] < -thre ) y[ind]<- eps else
              y[ind] <- exp(x[ind])
        }
      } else if ( "xtoa" == transform ) {
        for ( ind in seq_along(as.array(x)) ) {
          y[ind] <- .complexLog(x[ind])
        }
      } else if ( "gradfact" == transform )
        y <- x
    }
    else {
      limVal <- 36 
      y <- rep(0, length(x))
      switch (EXPR = transform,
              atox = { 
                index <- which(x< (-limVal)) 
                y[index] <- exp(-limVal) 
                x[index] <- NaN 
                index <- which(x<limVal) 
                y[index] <- exp(x[index]) 
                x[index] <- NaN 
                index <- which(!is.na(x)) 
                if (!(length(index) == 0))
                  y[index] <- exp(limVal) 
              }, 
              xtoa = {
                y <- log(x) }, 
              gradfact = {
                y <- x })
    }
    return (y)
  }
