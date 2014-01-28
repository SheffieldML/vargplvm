rbfard2VardistPsi0Compute <-
function (rbfardKern, vardist)
{
# % RBFARD2VARDISTPSI0COMPUTE  description not available.  
# % FORMAT
# % DESC 
# % description not available

k0 <- vardist$numData*rbfardKern$variance

return (k0)
}
