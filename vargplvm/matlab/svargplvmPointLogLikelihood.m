function [ll, model] = svargplvmPointLogLikelihood(model, vardistx, y)
% SVARGPLVMPOINTLOGLIKELIHOOD Log-likelihood of one or more points for the GP-LVM.
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
% ARG y : the NORMALISED observed data points for which the posterior is evaluated      
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
for i=model.testModalities
    indexMissing{i} = model.comp{i}.testPrecomp.indexMissing;
    indexPresent{i} = setdiff(1:model.comp{i}.d, indexMissing{i});
end
%y = y(:,indexPresent);   
    
if isfield(model, 'dynamics') && ~isempty(model.dynamics)
   dynUsed = 1;
else
   dynUsed = 0;
end
 
Nstar = size(y{1},1);
N = model.N;
% just keep the old model for possible future use 
%mOrig = model.m;  


% For the dynamics case we need to re-compute Psi statistics and other quantities.
% This is because of the coupling among the variational parameters; if one
% of the "free" means mu_bar change (or if new are added), then all original means need to be
% recalculated as mu=Kt*mu_bar. The same holds for covariances. And since
% the variational distribution changes, the Psi statistics and everything
% that depends on them needs to be re-computed.
if dynUsed == 1
    error('Not implemented for dynamics case')
end 

% LL1 TERM: This only takes into account the dimensions (m) from the training
% data (no test data at all here), where (m) is the set of dimensions
% which are missing in the test data.
%
ll1 = 0;
for i=model.testModalities
    if ~isempty(indexMissing)
        dmis = prod(size(indexMissing{i}));
        
        % Precompute again the parts that contain Y
        TrYY = model.comp{i}.testPrecomp.TrYY; %sum(sum(mOrig(:,indexMissing) .* mOrig(:,indexMissing)));
        
        ll1 = ll1 -0.5*(dmis*(-(model.N-model.comp{i}.k)*log(model.comp{i}.beta) ...
            + model.comp{i}.logDetAt) ...
            - (model.comp{i}.testPrecomp.TrPP ...
            - TrYY)*model.comp{i}.beta);
        ll1 = ll1 - 0.5*model.comp{i}.beta*dmis*model.comp{i}.Psi0 + 0.5*dmis*model.comp{i}.beta*model.comp{i}.TrC;
        ll1 = ll1-dmis*model.N/2*log(2*pi);
    end
end

% LL2 TERM : This takes into account the whole distribution of the training
% and test data but only for the dimensions which are partly observed in
% the test data. Now we need to calculate the parts of the Psi
% statistics that correspond to the test data and then amend the old Psi's
% (which were calculated above only for the training part) with these
% values.

model.N =  model.N + Nstar;

ll2 = 0;
for i=model.testModalities
    pointPsi0 = kernVardistPsi0Compute(model.comp{i}.kern, vardistx);
    pointPsi1 = kernVardistPsi1Compute(model.comp{i}.kern, vardistx, model.comp{i}.X_u);
    pointPsi2 = kernVardistPsi2Compute(model.comp{i}.kern, vardistx, model.comp{i}.X_u);
    
    model.comp{i}.N =  model.comp{i}.N + Nstar;
    model.comp{i}.d = prod(size(indexPresent{i}));
    model.comp{i}.Psi1 = [model.comp{i}.Psi1; pointPsi1];
    model.comp{i}.Psi2 = model.comp{i}.Psi2 + pointPsi2;
    model.comp{i}.Psi0 = model.comp{i}.Psi0 + pointPsi0;
    % normalize y exactly as model.m is normalized
    %my = y - repmat(model.bias(indexPresent),size(y,1),1);
    %my = my./repmat(model.scale(indexPresent),size(y,1),1);
    
    %
    % change the data (by including the new point and taking only the present indices)
    %model.m = mOrig(:,indexPresent);
    %model.m = [model.m; my];
    %model.TrYY = sum(sum(model.m .* model.m));
    model.comp{i}.TrYY = model.comp{i}.testPrecomp.TrYY2;
    model.comp{i}.C = model.comp{i}.invLm * model.comp{i}.Psi2 * model.comp{i}.invLmT;
    model.comp{i}.TrC = sum(diag(model.comp{i}.C)); % Tr(C)
    model.comp{i}.At = (1/model.comp{i}.beta) * eye(size(model.comp{i}.C,1)) + model.comp{i}.C; % At = beta^{-1} I + C
    model.comp{i}.Lat = jitChol(model.comp{i}.At)';
    model.comp{i}.invLat = model.comp{i}.Lat\eye(size(model.comp{i}.Lat,1));
    model.comp{i}.invLatT = model.comp{i}.invLat';
    model.comp{i}.logDetAt = 2*(sum(log(diag(model.comp{i}.Lat)))); % log |At|
    model.comp{i}.P1 = model.comp{i}.invLat * model.comp{i}.invLm; % M x M
    if isfield(model.comp{i}, 'DgtN') && model.comp{i}.DgtN
        model.comp{i}.P = model.comp{i}.P1 * (model.comp{i}.Psi1' * model.comp{i}.testPrecomp.mReduced);
    else
        % y is already normalised from svargplvmOptimisePoint
        model.comp{i}.P = model.comp{i}.P1 * (model.comp{i}.Psi1' *  [model.comp{i}.m(:, indexPresent{i}); y{i}]);
    end
    model.comp{i}.TrPP = sum(sum(model.comp{i}.P .* model.comp{i}.P));
    
    ll2 = ll2 -0.5*(model.comp{i}.d*(-(model.comp{i}.N-model.comp{i}.k)*log(model.comp{i}.beta) ...
        + model.comp{i}.logDetAt) ...
        - (model.comp{i}.TrPP ...
        - model.comp{i}.TrYY)*model.comp{i}.beta);
    ll2 = ll2 - 0.5*model.comp{i}.beta*model.comp{i}.d*model.comp{i}.Psi0 ...
        + 0.5*model.comp{i}.d*model.comp{i}.beta*model.comp{i}.TrC;
    
    ll2 = ll2-model.comp{i}.d*model.comp{i}.N/2*log(2*pi);
end

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

