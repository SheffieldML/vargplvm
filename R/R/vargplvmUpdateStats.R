vargplvmUpdateStats <-
function (model, X_u)
{
# % VARGPLVMUPDATESTATS Update stats for VARGPLVM model.
# % FORMAT
# % DESC no description available.
# % COPYRIGHT : Michalis K. Titsias, 2009-2011
# % COPYRIGHT : Neil D. Lawrence, 2009-2011
# % COPYRIGHT : Andreas C. Damianou, 2010-2011
# % 
# % SEEALSO : vargplvmOptimise, vargplvmExpandParam
# 
# % VARGPLVM
  
jitter <- 1e-6 
# %model.jitter = 1e-6 


# % %%% Precomputations for the KL term %%%
# %%% Maybe we should add something like model.dynamics.X =
# %%% model.dynamics.vardist.means (if visualisation requires that).
if (("dynamics" %in% names(model)) && (length(model$dynamics) > 0))
     model <- vargplvmDynamicsUpdateStats(model)  # not implemented

# %%% Precomputations for (the likelihood term of) the bound %%%

model$K_uu <- kernCompute(model$kern, X_u) 

# % Always add jitter (so that the inducing variables are "jitter" function variables)
# % and the above value represents the minimum jitter value
# % Putting jitter always ("if" in comments) is like having a second
# % whiteVariance in the kernel which is constant.
# %if (~isfield(model$kern, 'whiteVariance')) | model$kern$whiteVariance < jitter
#    % There is no white noise term so add some jitter.
if (!(model$kern$type == "rbfardjit"))
  model$K_uu <- model$K_uu + diag(rep(jitter, dim(model$K_uu)[1])) 


model$Psi0 <- kernVardistPsi0Compute(model$kern, model$vardist)
model$Psi1 <- (kernVardistPsi1Compute(model$kern, model$vardist, X_u))$Psi1
temp <- kernVardistPsi2Compute(model$kern, model$vardist, X_u)
model$Psi2 <- temp$Psi2
AS <- temp$P

# %K_uu_jit <- model$K_uu + model$jitter*eye(model$k) 
# %model$Lm <- chol(K_uu_jit, 'lower')           
  
#                                       % M is model$k
# %model$Lm <- chol(model$K_uu, 'lower')  

model$Lm <- t(jitChol(model$K_uu))       #% M x M: L_m (lower triangular)   ---- O(m^3)
model$invLm <- solve(model$Lm) %*% diag(model$k)    #% M x M: L_m^{-1}                 ---- O(m^3)

model$invLmT <- t(model$invLm)  #% L_m^{-T}
model$C <- (model$invLm %*% model$Psi2) %*% model$invLmT 

model$TrC <- sum(diag(model$C))  #% Tr(C)
# % Matrix At replaces the matrix A of the old implementation;  At is more stable
# % since it has a much smaller condition number than A=sigma^2 K_uu + Psi2
model$At <- diag((1/model$beta), dim(model$C)[1], dim(model$C)[1]) + model$C  #% At = beta^{-1} I + C

model$Lat <- t(jitChol(model$At)) 

model$invLat <- solve(model$Lat) %*% diag(dim(model$Lat)[1])

model$invLatT <- t(model$invLat) 
model$logDetAt <- 2*(sum(log(diag(model$Lat))))  #% log |At|

model$P1 <- model$invLat %*% model$invLm  #% M x M

# % First multiply the two last factors  so, the large N is only involved
# % once in the calculations (P1: MxM, Psi1':MxN, Y: NxD)
model$P <- model$P1 %*% (t(model$Psi1) %*% model$m) 

# % Needed for both, the bound's and the derivs. calculations.
model$TrPP <- sum(model$P * model$P)


# %%% Precomputations for the derivatives (of the likelihood term) of the bound %%%

# %model$B <- model$invLmT * model$invLatT * model$P  %next line is better
model$B <- t(model$P1) %*% model$P 
model$invK_uu <- model$invLmT %*% model$invLm 

Tb <- c((1/model$beta) * model$d) * (t(model$P1) %*% model$P1)
Tb <- Tb + (model$B %*% t(model$B)) 
model$T1 <- model$d * model$invK_uu - Tb 
model$X <- model$vardist$means 


return (model)
}
