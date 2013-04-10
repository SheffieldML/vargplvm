function [g, model] = vargplvmSeqDynLogLikeGradient(model, vardistx, y)
% VARGPLVMSEQDYNLOGLIKEGRADIENT Log-likelihood gradient for of a point of the GP-LVM.
% FORMAT
% DESC returns the gradient of the log likelihood with respect to
% the latent position, where the log likelihood is conditioned on
% the training set. 
% ARG model : the model for which the gradient computation is being
% done.
% ARG x : the latent position where the gradient is being computed.
% ARG y : the position in data space for which the computation is
% being done.
% RETURN g : the gradient of the log likelihood, conditioned on the
% training data, with respect to the latent position.
%
% SEEALSO : vargplvmPointLogLikelihood, vargplvmOptimisePoint, vagplvmSequenceLogLikeGradient
%
% COPYRIGHT : Michalis K. Titsias and Andreas Damianou, 2011

% VARGPLVM


% y is a new block/sequence 
Nstar = size(y,1);
N = model.N;

mask = sum(isnan(y),2); 
indexObservedData = find(mask==0)'; 
indexMissingData = setdiff(1:Nstar, indexObservedData);

% Compute fully observed test data points and partially 
% observed data points 
yOb = y(indexObservedData, :);
yMs = y(indexMissingData, :);

% Indices of missing dimension in the Missingdata
indexMissing = [];
indexPresent = [1:model.d];
if ~isempty(yMs)
   indexMissing = find(isnan(yMs(1,:)));
   indexPresent = setdiff(1:model.d, indexMissing);
   yMs = yMs(:,indexPresent);   
end
    
% normalize yOb and yMs exactly as model.m is normalized 
myOb = yOb;
if ~isempty(yOb)
  myOb = yOb - repmat(model.bias,size(yOb,1),1); 
  myOb = myOb./repmat(model.scale,size(yOb,1),1); 
end
myMs = yMs;
if ~isempty(yMs)
   myMs = yMs - repmat(model.bias(indexPresent),size(yMs,1),1);  
   myMs = myMs./repmat(model.scale(indexPresent),size(yMs,1),1);  
end
mOrig = model.m;
%modelOrig = model;

% re-order test data so that observed are first and then are the missing 
Order = [indexObservedData, indexMissingData];
model.dynamics.t_star = model.dynamics.t_star(Order);
vardistx.means = vardistx.means(Order,:);
vardistx.covars = vardistx.covars(Order,:);
nObsData = size(indexObservedData,2);

% Form the modelTest for the new block which allows to compute the variational
% distribution (given the barmus and lambdas) for this new block. Notice that this distribution
% is independent from the other blocks in the training set. 
modelTest = model;
%modelTest.m = [myOb; myMs]; 
modelTest.dynamics.t = model.dynamics.t_star;
modelTest.dynamics.vardist = vardistx; 
modelTest.dynamics.N = size(modelTest.dynamics.vardist.means, 1);
modelTest.vardist.numData = modelTest.dynamics.N;
modelTest.vardist.nParams = 2*prod(size(modelTest.dynamics.vardist.means));
Kt = kernCompute(model.dynamics.kern, model.dynamics.t_star);
modelTest.dynamics.Kt = Kt;
modelTest.dynamics.fixedKt = 1;
% This enables the computation of the variational means and covars
modelTest.dynamics.seq = []; 
modelTest = vargpTimeDynamicsUpdateStats(modelTest);


% The fully observed subset of data from the new block must be augmented with all previous
% training data to form the LL1 term in the whole bound 
gVarmeansLik = zeros(1, nObsData*model.dynamics.q);
gVarcovsLik = zeros(1, nObsData*model.dynamics.q);
if ~isempty(yOb)
   vardistOb = modelTest.vardist; 
   vardistOb.means = modelTest.vardist.means(1:nObsData, :);
   vardistOb.covars = modelTest.vardist.covars(1:nObsData, :);
   vardistOb.nParams = 2*prod(size(vardistOb.means));
   vardistOb.numData = size(vardistOb.means,1);
 
   % Psi statistics for the data of the new block/sequence which have a fully
   % observed features/dimensions
   obsPsi0 = kernVardistPsi0Compute(model.kern, vardistOb);
   obsPsi1 = kernVardistPsi1Compute(model.kern, vardistOb, model.X_u);
   obsPsi2 = kernVardistPsi2Compute(model.kern, vardistOb, model.X_u);
   
   model.N = N + size(yOb,1);
   model.m = [model.m; myOb];
   model.d = prod(size(indexMissing));
   model.Psi1 = [model.Psi1; obsPsi1]; 
   model.Psi2 = model.Psi2 + obsPsi2;
   model.Psi0 = model.Psi0 + obsPsi0; 
   model.C = model.invLm * model.Psi2 * model.invLmT;
   model.At = (1/model.beta) * eye(size(model.C,1)) + model.C; % At = beta^{-1} I + C
   model.Lat = jitChol(model.At)';
   model.invLat = model.Lat\eye(size(model.Lat,1));  
   model.P1 = model.invLat * model.invLm; % M x M
   P1TP1 = (model.P1' * model.P1);
   model.P = model.P1 * (model.Psi1' * model.m(:,indexMissing));
   model.B = model.P1' * model.P;
   Tb = (1/model.beta) * model.d * (model.P1' * model.P1);
        Tb = Tb + (model.B * model.B');
   model.T1 = model.d * model.invK_uu - Tb;
  
   % Precompuations for the gradients
   gPsi1 = model.beta*(P1TP1*model.Psi1'*model.m(:,indexMissing)*myOb(:,indexMissing)');
   %gPsi1 = model.beta * model.m(:,indexMissing) * model.B';
   %gPsi1 = gPsi1'; % because it is passed to "kernVardistPsi1Gradient" as gPsi1'...
   gPsi2 = (model.beta/2) * model.T1;
   gPsi0 = -0.5 * model.beta * model.d;
   [gKern1, gVarmeans1, gVarcovs1, gInd1] = kernVardistPsi1Gradient(model.kern, vardistOb, model.X_u, gPsi1');
   [gKern2, gVarmeans2, gVarcovs2, gInd2] = kernVardistPsi2Gradient(model.kern, vardistOb, model.X_u, gPsi2);
   [gKern0, gVarmeans0, gVarcovs0] = kernVardistPsi0Gradient(model.kern, vardistOb, gPsi0);
   gVarmeansLik = gVarmeans0 + gVarmeans1 + gVarmeans2;
   gVarcovsLik = gVarcovs0 + gVarcovs1 + gVarcovs2;   
%
end 


% GRADIENT FOR THE LL1 TERM
%
gPointDyn1 = zeros(Nstar, model.dynamics.q*2);
if ~isempty(indexMissing)
%    
   % means
   gPointDyn1(:,1:model.dynamics.q) = modelTest.dynamics.Kt(1:nObsData,:)'*reshape(gVarmeansLik, nObsData, model.q);
   % covars
   gVarcovsLik = reshape(gVarcovsLik, nObsData, model.q);
   gcovTmp = zeros(Nstar,model.dynamics.q);
   for q=1:model.dynamics.q
       LambdaH_q = modelTest.dynamics.vardist.covars(:,q).^0.5;
       Bt_q = eye(Nstar) + LambdaH_q*LambdaH_q'.*modelTest.dynamics.Kt;
       % Invert Bt_q
       Lbt_q = jitChol(Bt_q)';
       G1 = Lbt_q \ diag(LambdaH_q);
       G = G1*modelTest.dynamics.Kt;
       % Find Sq
       Sq = modelTest.dynamics.Kt - G'*G;
       Sq = - (Sq .* Sq);
       % only the cross matrix is needed as in barmu case
       gcovTmp(:,q) = Sq(1:nObsData,:)'*gVarcovsLik(:,q);
   end
   gPointDyn1(:,(model.dynamics.q+1):(model.dynamics.q*2)) = gcovTmp.*modelTest.dynamics.vardist.covars;
%   
end


% GRADIENT FOR THE LL2 TERM PLUS KL TERM
%
% The data in the test block/seq with missing values are processed to get their
% psi statisitcs. This is needed to form the LL2 term that is the part 
% of the bound corresponding to dimensions observed everywhere (in all
% training and test points)
vardistMs = modelTest.vardist; 
vardistMs.means = modelTest.vardist.means(nObsData+1:end, :);
vardistMs.covars = modelTest.vardist.covars(nObsData+1:end, :);
vardistMs.nParams = 2*prod(size(vardistMs.means));
vardistMs.numData = size(vardistMs.means,1);
% Psi statistics for the data of the new block/sequence which have partially
% observed dimensions
missingPsi0 = kernVardistPsi0Compute(model.kern, vardistMs);
missingPsi1 = kernVardistPsi1Compute(model.kern, vardistMs, model.X_u);
missingPsi2 = kernVardistPsi2Compute(model.kern, vardistMs, model.X_u);  
model.m = [mOrig(:,indexPresent); myOb(:, indexPresent); myMs];  
model.N = N + Nstar;
model.d = prod(size(indexPresent));
model.dynamics.nParams = model.dynamics.nParams + 2*prod(size(vardistx.means));
model.nParams = model.nParams + 2*prod(size(vardistx.means));
model.Psi1 = [model.Psi1; missingPsi1]; 
model.Psi2 = model.Psi2 + missingPsi2;
model.Psi0 = model.Psi0 + missingPsi0;
model.C = model.invLm * model.Psi2 * model.invLmT;
model.At = (1/model.beta) * eye(size(model.C,1)) + model.C; % At = beta^{-1} I + C
model.Lat = jitChol(model.At)';
model.invLat = model.Lat\eye(size(model.Lat,1));  
model.P1 = model.invLat * model.invLm; % M x M
P1TP1 = (model.P1' * model.P1);
model.P = model.P1 * (model.Psi1' * model.m);
model.B = model.P1' * model.P;
Tb = (1/model.beta) * model.d * (model.P1' * model.P1);
     Tb = Tb + (model.B * model.B');
model.T1 = model.d * model.invK_uu - Tb;
     
% Precompuations for the gradients
gPsi1 = model.beta*(P1TP1*model.Psi1'*model.m*[myOb(:,indexPresent); myMs]');
gPsi2 = (model.beta/2) * model.T1;
gPsi0 = -0.5 * model.beta * model.d;
[gKern1, gVarmeans1, gVarcovs1, gInd1] = kernVardistPsi1Gradient(model.kern, modelTest.vardist, model.X_u, gPsi1');
[gKern2, gVarmeans2, gVarcovs2, gInd2] = kernVardistPsi2Gradient(model.kern, modelTest.vardist, model.X_u, gPsi2);
[gKern0, gVarmeans0, gVarcovs0] = kernVardistPsi0Gradient(model.kern, modelTest.vardist, gPsi0);
gVarmeansLik = gVarmeans0 + gVarmeans1 + gVarmeans2;
gVarcovsLik = gVarcovs0 + gVarcovs1 + gVarcovs2;      
%

gPointDyn2 = zeros(Nstar, model.dynamics.q*2); 
% means
gPointDyn2(:,1:model.dynamics.q) = modelTest.dynamics.Kt'*(reshape(gVarmeansLik, Nstar, model.q) - modelTest.dynamics.vardist.means); 
% covars
gVarcovsLik = reshape(gVarcovsLik, Nstar, model.q);
gcovTmp = zeros(Nstar,model.dynamics.q);
for q=1:model.dynamics.q
    LambdaH_q = modelTest.dynamics.vardist.covars(:,q).^0.5;
    Bt_q = eye(Nstar) + LambdaH_q*LambdaH_q'.*modelTest.dynamics.Kt;
    % Invert Bt_q
    Lbt_q = jitChol(Bt_q)';
    G1 = Lbt_q \ diag(LambdaH_q);
    G = G1*modelTest.dynamics.Kt;
    % Find Sq
    Sq = modelTest.dynamics.Kt - G'*G;
    gcovTmp(:,q) = - (Sq .* Sq) * (gVarcovsLik(:,q) + 0.5*modelTest.dynamics.vardist.covars(:,q));
end
gPointDyn2(:,(model.dynamics.q+1):(model.dynamics.q*2)) = gcovTmp.*modelTest.dynamics.vardist.covars; 

gPointDyn = gPointDyn1 + gPointDyn2;

% applying the inverse Order to give the gradeitn in the 
%orginal order of the data
T = eye(Nstar);
T = T(Order,:);
% inverse permutation
T = T';
% there must be a better way to do this
InvOrder = zeros(1,Nstar);
for i=1:Nstar
    InvOrder(i) = find(T(i,:)==1); 
end
gPointDyn = gPointDyn(InvOrder, :);

g = gPointDyn(:)';
