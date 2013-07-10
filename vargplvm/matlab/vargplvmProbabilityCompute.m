function logProb = vargplvmProbabilityCompute(model, y, display, iters)

% VARGPLVMPROBABILITYCOMPUTE description
  
% VARGPLVM
  
% Takes an input a trained vargplvm and a test data point (with possibly missing values)
% Computes the probability density in the test data point 


% Indices of missing dimension
indexMissing = find(isnan(y(1,:)));
indexPresent = setdiff(1:model.d, indexMissing);
y = y(:,indexPresent); 


% compute the variational lower without the new data point 
Fold = vargplvmLogLikelihood(model);

% initialize the latent point using the nearest neighbour 
% from the training data
dst = dist2(y(indexPresent), model.y(:,indexPresent));
[mind, mini] = min(dst);
% create the variational distribtion for the test latent point
vardistx = vardistCreate(model.vardist.means(mini,:), model.q, 'gaussian');

% optimize over the latent point 
model.vardistx = vardistx;
[X, varX] = vargplvmOptimisePoint(model, vardistx, y, display, iters);

% compute the variational lower with the new data point included
vardistx = model.vardistx;
vardistx.means = X; 
vardistx.covars = varX;
Fnew = vargplvmPointLogLikelihood(model, vardistx, y);

% compute the probability 
logProb = Fnew - Fold;
 