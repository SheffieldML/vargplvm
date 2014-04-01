rbfard2biasVardistPsi2Compute <-
function(rbfardKern, biasKern, vardist, Z)
{
# % RBFARD2BIASVARDISTPSI2COMPUTE  description not available.  
# % FORMAT
# % DESC 
# % description not available

Psi1 <- (rbfard2VardistPsi1Compute(rbfardKern, vardist, Z))$K 

sumPsi1 <- colSums(Psi1)  

Psi2 <- matrix(1, dim(Z)[1],1)%*%sumPsi1  

Pnobias <- Psi2 + t(Psi2)  
Psi2 <- biasKern$variance*Pnobias  

# % second naive way
# %Psi22 <- zeros(size(Z,1),size(Z,1)) 
# %for j<-1:size(Z,1)
# %    for i<-1:size(Z,1)
# %        Psi22(j,i) <- biasKern.variance*(sum(Psi1(:,j)) + sum(Psi1(:,i))) 
# %    end
# %end
# %sum(sum(abs(Psi2 - Psi22))) 

return (list(Psi2 = Psi2, Pnobias = Pnobias, Psi1 = Psi1))
}
