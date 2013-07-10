function [mu, varsigma] = vargplvmPosteriorMeanVar(model, X, varX)

% VARGPLVMPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.
% FORMAT
% DESC returns the posterior mean and variance for a given set of
% points.
% ARG model : the model for which the posterior will be computed.
% ARG X : variational mean in the latent space for which posterior is computed.
% ARG varX : variational variances in the latent space for which posterior is computed (assumed zero if not present).
% RETURN mu : the mean of the posterior distribution.
% RETURN sigma : the variances of the posterior distributions.
%
% SEEALSO : gpPosteriorMeanVar, vargplvmCreate
%
% COPYRIGHT : Michalis K. Titsias and Neil D. Lawrence, 2009, 2011

% VARGPLVM


% do prediction by replacing the variational distribution with a delta function  
%model.K_uf = kernCompute(model.kern, model.X_u, model.vardist.means);
%model.A = (1/model.beta)*model.K_uu + model.K_uf*model.K_uf';
%[model.Ainv, U] = pdinv(model.A);
%[mu1, varsigma1] = gpPosteriorMeanVar(model, vardistX.means);


% Find exactly the mean and the variances of the predictive distribution
% (which is not Gaussian, however its moments can be computed in closed-form)

if nargin < 3
  vardistX.covars = repmat(0.0, size(X, 1), size(X, 2));%zeros(size(X, 1), size(X, 2));
else
  vardistX.covars = varX;
end
vardistX.latentDimension = size(X, 2);
vardistX.numData = size(X, 1);
%model.vardist.covars = 0*model.vardist.covars; 
vardistX.means = X;
%model = vargplvmUpdateStats(model, model.X_u);


Ainv = model.P1' * model.P1; % size: NxN

if ~isfield(model,'alpha')
    model.alpha = Ainv*model.Psi1'*model.m; % size: 1xD
end
Psi1_star = kernVardistPsi1Compute(model.kern, vardistX, model.X_u);

% mean prediction 
mu = Psi1_star*model.alpha; % size: 1xD

if nargout > 1
   % 
   % precomputations
   vard = vardistCreate(zeros(1,model.q), model.q, 'gaussian');
   Kinvk = (model.invK_uu - (1/model.beta)*Ainv);
   %
   for i=1:size(vardistX.means,1)
      %
      vard.means = vardistX.means(i,:);
      vard.covars = vardistX.covars(i,:);
      % compute psi0 term
      Psi0_star = kernVardistPsi0Compute(model.kern, vard);
      % compute psi2 term
      Psi2_star = kernVardistPsi2Compute(model.kern, vard, model.X_u);
    
      vars = Psi0_star - sum(sum(Kinvk.*Psi2_star));
      
      for j=1:model.d
         %[model.alpha(:,j)'*(Psi2_star*model.alpha(:,j)), mu(i,j)^2]
         varsigma(i,j) = model.alpha(:,j)'*(Psi2_star*model.alpha(:,j)) - mu(i,j)^2;  
      end
      varsigma(i,:) = varsigma(i,:) + vars; 
      %
   end
   % 
   if isfield(model, 'beta')
      varsigma = varsigma + (1/model.beta);
   end
   %
end
      
% Rescale the mean
mu = mu.*repmat(model.scale, size(vardistX.means,1), 1);

% Add the bias back in
mu = mu + repmat(model.bias, size(vardistX.means,1), 1);

% rescale the variances
if nargout > 1
    varsigma = varsigma.*repmat(model.scale.*model.scale, size(vardistX.means,1), 1);
end
  