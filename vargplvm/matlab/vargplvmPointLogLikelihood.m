function [ll, model] = vargplvmPointLogLikelihood(model, vardistx, y)
% VARGPLVMPOINTLOGLIKELIHOOD Log-likelihood of one or more points for the GP-LVM.
% FORMAT
% DESC returns the log likelihood of some latent points and the corresponding
% (possibly partly) observed data point for the posterior prediction of the GP-LVM model.
% It is perfectly functional when the input is just for one point. Note
% that this function assumes that vargplvmUpdateStats is not called before
% (unlike the vargplvmLogLikelihood function), so all the necessary
% precomputations (due to the change in the parameters during optimisation)
% are being done here. 
% ARG model : the model for which the points prediction will be
% made.
% ARG vardistx : the variational distribution over latent points for which the posterior distribution
% will be evaluated. It contains the mean and tha diagonal covarriance 
% ARG y : the observed data points for which the posterior is evaluated      
%
% SEEALSO : vargplvmCreate, vargplvmOptimisePoint, vargplvmPointObjective
%
% COPYRIGHT : Michalis K. Titsias and Neil D. Lawrence, 2009, 2011
% COPYRIGHT : Andreas C. Damianou, 2011
% VARGPLVM


% !!!!!! this function can become faster with precomputations stored in the
% structure model !!!!! 


% compute firstly the lower bound that corresponds to the missing indices
% in y 
%  --- this bit of the lower bound does not depend on the vardistx    

jitter = 1e-6;


% Indices of missing dimension
indexMissing = find(isnan(y(1,:)));
indexPresent = setdiff(1:model.d, indexMissing);
y = y(:,indexPresent);   
    
if isfield(model, 'dynamics') && ~isempty(model.dynamics)
   dynUsed = 1;
else
   dynUsed = 0;
end
 
Nstar = size(y,1);
N = model.N;
% just keep the old model for possible future use 
mOrig = model.m;  

Kt = zeros(N+Nstar, N+Nstar);

% For the dynamics case we need to re-compute Psi statistics and other quantities.
% This is because of the coupling among the variational parameters; if one
% of the "free" means mu_bar change (or if new are added), then all original means need to be
% recalculated as mu=Kt*mu_bar. The same holds for covariances. And since
% the variational distribution changes, the Psi statistics and everything
% that depends on them needs to be re-computed.
if dynUsed == 1
    %%% RE-OPT-CODE-NEW_
    if isfield(model.dynamics, 'reoptimise') && model.dynamics.reoptimise   
        % In this case the inducing points and theta_t are optimised and 
        % since updateStats is not called, any computations which need X_u
        % or theta_t must be done here.
        model = vargplvmDynamicsUpdateStats(model); % Update Kt, mu and S
        model.K_uu = kernCompute(model.kern, model.X_u); 
        model.K_uu = model.K_uu ...
            + sparseDiag(repmat(jitter, size(model.K_uu, 1), 1));
        model.Lm = jitChol(model.K_uu)';
        model.invLm = model.Lm\eye(model.k);
        model.invLmT = model.invLm';
    end %%% _RE-OPT-CODE-NEW 
    
   %model.y = model.y(:,indexMissing);
   model.m = model.m(:,indexMissing);
   model.d = prod(size(indexMissing));
   % Augment the time vector to include the timestamp of the new point
   model.dynamics.t = [model.dynamics.t; model.dynamics.t_star];
   % Augment the reparametrized variational parameters mubar and lambda
   model.dynamics.vardist.means = [model.dynamics.vardist.means; vardistx.means];
   model.dynamics.vardist.covars = [model.dynamics.vardist.covars; vardistx.covars]; 
   model.dynamics.N =  model.dynamics.N + Nstar;
   model.vardist.numData = model.dynamics.N;
   model.dynamics.vardist.numData = model.dynamics.N;
   model.vardist.nParams = 2*prod(size(model.dynamics.vardist.means));
   model.dynamics.vardist.nParams = 2*prod(size(model.dynamics.vardist.means));
   %model.dynamics.seq = model.dynamics.N; %%%%%% Chec
  
   % Faster computation of the Kt matrix 
   % (all this could have been avoided as Kt is constant... we need to do these
   % precomputations outside of this function) 
   % construct Kt 
   Kt(1:N,1:N) = model.dynamics.Kt;
   % new diagonal block
   Kt(N+1:end, N+1:end) = kernCompute(model.dynamics.kern, model.dynamics.t_star);
   % cross block; it only has to be computed if the model is not learning
   % from individual sequences (in which case, no correlations between
   % blocks need to be captured).
   if ~isfield(model.dynamics, 'seq') || isempty(model.dynamics.seq)
        Kt(1:N, N+1:end) = kernCompute(model.dynamics.kern, model.dynamics.t(1:N), model.dynamics.t_star);
        Kt(N+1:end, 1:N) = Kt(1:N, N+1:end)';
   end
   model.dynamics.Kt = Kt;
   % The following means that Kt needs not be recomputed from the input
   % vector t, because it is computed externally using the old Kt and some
   % additional columns/rows which augment it.
   model.dynamics.fixedKt = 1;
   % The following updates the variational distribution with the true
   % parameters and (if model.dynamics.fixedKt==0) recomputes Kt.
   model = vargpTimeDynamicsUpdateStats(model);
   %model.dynamics.fixedKt = 0;
   %model = vargpTimeDynamicsUpdateStats(model);
   
   % Temporary variational structure which holds the recomputed var.
   % parameters but only for the part that corresponds to the training
   % data. But since this depends on the test data as well (because of the
   % coupling) we have to calculate everything together (as we did above)
   % and then select only the rows that correspond to the tr. data.
   vardist2 = model.vardist;
   vardist2.means = model.vardist.means(1:end-Nstar,:);
   vardist2.covars = model.vardist.covars(1:end-Nstar,:);
   vardist2.nParams = 2*prod(size(vardist2.means));
   vardist2.numData = size(vardist2.means,1);
   model.Psi0 = kernVardistPsi0Compute(model.kern, vardist2);
   model.Psi1 = kernVardistPsi1Compute(model.kern, vardist2, model.X_u);
   model.Psi2 = kernVardistPsi2Compute(model.kern, vardist2, model.X_u);
    
   model.C = model.invLm * model.Psi2 * model.invLmT;
   model.TrC = sum(diag(model.C)); % Tr(C)
   model.At = (1/model.beta) * eye(size(model.C,1)) + model.C; % At = beta^{-1} I + C 
   model.Lat = jitChol(model.At)';
   model.invLat = model.Lat\eye(size(model.Lat,1));  
   model.invLatT = model.invLat';
   model.logDetAt = 2*(sum(log(diag(model.Lat)))); % log |At|
   model.P1 = model.invLat * model.invLm; % M x M
   model.P = model.P1 * (model.Psi1' * model.m);
   TrPP = sum(sum(model.P .* model.P));
else
   P = model.P1 * (model.Psi1' * model.m(:,indexMissing));
   TrPP = sum(sum(P .* P));
end 



% LL1 TERM: This only takes into account the dimensions (m) from the training
% data (no test data at all here), where (m) is the set of dimensions
% which are missing in the test data.
%
ll1 = 0;
if ~isempty(indexMissing)
   dmis = prod(size(indexMissing));
   
   % Precompute again the parts that contain Y
   TrYY = sum(sum(mOrig(:,indexMissing) .* mOrig(:,indexMissing)));
   
   ll1 = -0.5*(dmis*(-(model.N-model.k)*log(model.beta) ...
				  + model.logDetAt) ...
	      - (TrPP ...
	      - TrYY)*model.beta);
   ll1 = ll1 - 0.5*model.beta*dmis*model.Psi0 + 0.5*dmis*model.beta*model.TrC;
   ll1 = ll1-dmis*model.N/2*log(2*pi);
end

% LL2 TERM : This takes into account the whole distribution of the training
% and test data but only for the dimensions which are partly observed in
% the test data. Now we need to calculate the parts of the Psi
% statistics that correspond to the test data and then amend the old Psi's
% (which were calculated above only for the training part) with these
% values.

if dynUsed  
    vardist2 = model.vardist;
    vardist2.means = model.vardist.means(end-Nstar+1:end,:);
    vardist2.covars = model.vardist.covars(end-Nstar+1:end,:);
    vardist2.nParams = 2*prod(size(vardist2.means));
    vardist2.numData = size(vardist2.means,1);
    pointPsi0 = kernVardistPsi0Compute(model.kern, vardist2);
    pointPsi1 = kernVardistPsi1Compute(model.kern, vardist2, model.X_u);
    pointPsi2 = kernVardistPsi2Compute(model.kern, vardist2, model.X_u);
else
    pointPsi0 = kernVardistPsi0Compute(model.kern, vardistx);
    pointPsi1 = kernVardistPsi1Compute(model.kern, vardistx, model.X_u);
    pointPsi2 = kernVardistPsi2Compute(model.kern, vardistx, model.X_u);
end

model.N =  model.N + Nstar;
model.d = prod(size(indexPresent));
model.Psi1 = [model.Psi1; pointPsi1]; 
model.Psi2 = model.Psi2 + pointPsi2;
model.Psi0 = model.Psi0 + pointPsi0;
% normalize y exactly as model.m is normalized 
my = y - repmat(model.bias(indexPresent),size(y,1),1);
my = my./repmat(model.scale(indexPresent),size(y,1),1);



%
% change the data (by including the new point and taking only the present indices)
model.m = mOrig(:,indexPresent);
model.m = [model.m; my]; 
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

ll2 = -0.5*(model.d*(-(model.N-model.k)*log(model.beta) ...
				  + model.logDetAt) ...
	      - (model.TrPP ...
	      - model.TrYY)*model.beta);
ll2 = ll2 - 0.5*model.beta*model.d*model.Psi0 + 0.5*model.d*model.beta*model.TrC;

ll2 = ll2-model.d*model.N/2*log(2*pi);


% KL TERM 
%  
if dynUsed
    KLdiv = modelVarPriorBound(model);
else
    model.vardist.means = [model.vardist.means; vardistx.means];
    model.vardist.covars = [model.vardist.covars; vardistx.covars];
    varmeans = sum(sum(model.vardist.means.*model.vardist.means)); 
    varcovs = sum(sum(model.vardist.covars - log(model.vardist.covars)));
    KLdiv = -0.5*(varmeans + varcovs) + 0.5*model.q*model.N; 
end

% sum all terms
ll = ll1 + ll2 + KLdiv; 

