spdiags <-
function (arg1,arg2,arg3,arg4)
{
# % SPDIAGS description not available.  
# % FORMAT
# % DESC 
# % description not available

B <- arg1 
if (is.matrix(arg2))
d <- matrix(arg2,dim(arg2)[1]*dim(arg2)[2],1)
else
d <- arg2

p <- length(d) 

A <- sparseMatrix(i = 1:arg3, j = 1:arg3, x = 0, dims= c(arg3,arg4))
cat("in spdiags.R ")
print(A)
cat(" arg3 ")
print(arg3)
cat(" arg4 ")
pring(arg4)

m <- dim(A)[1] 
n <- dim(A)[2] 


len<-matrix(0,p+1,1)
for (k in 1:p)
len[k+1] <- len[k]+length(max(1,1-d[k]):min(m,n-d[k])) 

a <- matrix(0, len[p+1],3) 
for (k in 1:p)
{
# % Append new d(k)-th diagonal to compact form
i <- t(max(1,1-d[k]):min(m,n-d[k])) 
a[(len[k]+1):len[k+1],] <- c(i, i+d[k], B[(i+(m>=n)*d[k]),k]) 
}

res1 <- sparseMatrix(i = a[,1], j = a[,2], x = a[,3], dims = c(m,n)) 


return (res1)
}
