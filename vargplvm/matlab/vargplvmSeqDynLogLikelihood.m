function [ll, model] = vargplvmSeqDynLogLikelihood(model, vardistx, y)
% VARGPLVMSEQDYNLOGLIKELIHOOD Log-likelihood of a point for the GP-LVM.
% FORMAT
% DESC returns the log likelihood of a latent point and an observed
% data point for the posterior prediction of the GP-LVM model.
% ARG model : the model for which the point prediction will be
% made.
% ARG vardistx : the variational distribution over latent point for which the posterior distribution
% will be evaluated. It contains the mean and the diagonal covarriance 
% ARG y : the observed data point for which the posterior is evaluated
% ARG indexPresent: indicates which indices from the observed vector are present
%         (e.g. when when all indices are present, then y will d-dimensional     
%          and indexPresent = 1:D)          
%
% SEEALSO : vargplvmCreate, vargplvmOptimisePoint, vargplvmPointObjective
%
% COPYRIGHT : Michalis K. Titsias and Andreas Damianou, 2011

% VARGPLVM


% !!!!!! this function can become faster with precomputations stored in the
% structure model !!!!! 


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
modelTest.dynamics.N =  size(modelTest.dynamics.vardist.means, 1);
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
   
   model.N = model.N + size(yOb,1);
   model.m = [model.m; myOb];
   model.Psi1 = [model.Psi1; obsPsi1]; 
   model.Psi2 = model.Psi2 + obsPsi2;
   model.Psi0 = model.Psi0 + obsPsi0;
   model.C = model.invLm * model.Psi2* model.invLmT;
   model.TrC = sum(diag(model.C)); % Tr(C)
   model.At = (1/model.beta) * eye(size(model.C,1)) + model.C; % At = beta^{-1} I + C 
   model.Lat = jitChol(model.At)';
   model.invLat = model.Lat\eye(size(model.Lat,1));  
   model.invLatT = model.invLat';
   model.logDetAt = 2*(sum(log(diag(model.Lat)))); % log |At|
   model.P1 = model.invLat * model.invLm; % M x M
   model.P = model.P1 * (model.Psi1' * model.m(:,indexMissing));
   TrPP = sum(sum(model.P .* model.P));
else
   P = model.P1 * (model.Psi1' * model.m(:,indexMissing));
   TrPP = sum(sum(P .* P));
end 



% LL1 TERM
%
ll1 = 0;
if ~isempty(indexMissing)
   dmis = prod(size(indexMissing));
   
   % Precompute again the parts that contain Y
   TrYY = sum(sum(model.m(:,indexMissing) .* model.m(:,indexMissing)));
   
   ll1 = -0.5*(dmis*(-(model.N-model.k)*log(model.beta) ...
				  + model.logDetAt) ...
	      - (TrPP ...
	      - TrYY)*model.beta);
   ll1 = ll1 - 0.5*model.beta*dmis*model.Psi0 + 0.5*dmis*model.beta*model.TrC;
   ll1 = ll1-dmis*model.N/2*log(2*pi);
end


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
model.N =  N + Nstar;
model.d = prod(size(indexPresent));
model.Psi1 = [model.Psi1; missingPsi1]; 
model.Psi2 = model.Psi2 + missingPsi2;
model.Psi0 = model.Psi0 + missingPsi0;
model.TrYY = sum(sum(model.m .* model.m));
model.C = model.invLm * model.Psi2 * model.invLmT;
model.TrC = sum(diag(model.C)); % Tr(C)
model.At = (1/model.beta) * eye(size(model.C,1)) + model.C; % At = beta^{-1} I + C
model.Lat = jitChol(model.At)';
model.invLat = model.Lat\eye(size(model.Lat,1));  
model.invLatT = model.invLat';
model.logDetAt = 2*(sum(log(diag(model.Lat)))); % log |At|
model.P1 = model.invLat * model.invLm; % M x M
model.P = model.P1 * (model.Psi1' * model.m);
model.TrPP = sum(sum(model.P .* model.P));
    

% LL2 TERM 
%
ll2 = 0;
if ~isempty(indexPresent)
%   
    ll2 = -0.5*(model.d*(-(model.N-model.k)*log(model.beta) ...
				  + model.logDetAt) ...
	      - (model.TrPP ...
	      - model.TrYY)*model.beta);
    ll2 = ll2 - 0.5*model.beta*model.d*model.Psi0 + 0.5*model.d*model.beta*model.TrC;

    ll2 = ll2-model.d*model.N/2*log(2*pi);
%    
end


% KL TERM 
%
model.dynamics.t = [model.dynamics.t; model.dynamics.t_star];
% Augment the reparametrized variational parameters mubar and lambda
model.dynamics.vardist.means = [model.dynamics.vardist.means; vardistx.means];
model.dynamics.vardist.covars = [model.dynamics.vardist.covars; vardistx.covars]; 
model.dynamics.N = N + Nstar;
model.vardist.numData = model.dynamics.N;
Kt = zeros(N+Nstar,N+Nstar);
Kt(1:N,1:N) = model.dynamics.Kt; 
Kt(N+1:end, N+1:end) = modelTest.dynamics.Kt; 
model.dynamics.Kt = Kt;
model.dynamics.fixedKt = 1;
model.dynamics.seq = [model.dynamics.seq, (N+Nstar)];
KLdiv = modelVarPriorBound(model);

% sum all terms
ll = ll1 + ll2 + KLdiv; 
