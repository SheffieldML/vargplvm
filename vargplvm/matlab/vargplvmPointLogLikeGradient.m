function [g, model] = vargplvmPointLogLikeGradient(model, vardistx, y)
% VARGPLVMPOINTLOGLIKEGRADIENT Log-likelihood gradient for of some points of the GP-LVM.
% FORMAT
% DESC returns the gradient of the log likelihood with respect to
% the latent positions, where the log likelihood is conditioned on
% the training set. The function works even when we are interested only in
% one point. See vargplvmPointLogLikelihood for more details.
% ARG model : the model for which the gradient computation is being
% done.
% ARG x : the latent positions where the gradient is being computed.
% ARG y : the positions in data space for which the computation is
% being done.
% RETURN g : the gradient of the log likelihood, conditioned on the
% training data, with respect to the latent position.
%
% SEEALSO : vargplvmPointLogLikelihood, vargplvmOptimisePoint, vagplvmSequenceLogLikeGradient
%
% COPYRIGHT : Michalis K. Titsias and Neil D. Lawrence, 2009, 2011
% COPYRIGHT : Andreas C. Damianou, 2011
% VARGPLVM


if isfield(model, 'dynamics') && ~isempty(model.dynamics)
    dynUsed = 1;
else
    dynUsed = 0;
end

indexMissing = find(isnan(y(1,:)));
indexPresent = setdiff(1:model.d, indexMissing);
y = y(:,indexPresent);


if dynUsed
    % use a special function
    [g, model] = dynPointLogLikeGradient(model, vardistx, y, indexPresent, indexMissing);
    % g = g(:)'; %%% RE-OPT-CODE-REM
    return;
end



% normalize y exactly as model.m is normalized
my = y - repmat(model.bias(indexPresent),size(y,1),1);
my = my./repmat(model.scale(indexPresent),size(y,1),1);

[gPsi0, gPsi1, gPsi2] = vargpCovGrads(model, vardistx, my, indexPresent);

[gKern1, gVarmeans1, gVarcovs1] = kernVardistPsi1Gradient(model.kern, vardistx, model.X_u, gPsi1');
[gKern2, gVarmeans2, gVarcovs2] = kernVardistPsi2Gradient(model.kern, vardistx, model.X_u, gPsi2);
[gKern0, gVarmeans0, gVarcovs0] = kernVardistPsi0Gradient(model.kern, vardistx, gPsi0);

gVarcovs0 = (gVarcovs0(:).*vardistx.covars(:))';
gVarcovs1 = (gVarcovs1(:).*vardistx.covars(:))';
gVarcovs2 = (gVarcovs2(:).*vardistx.covars(:))';

% KL divergence terms
gVarmeansKL = - vardistx.means(:)';
% !!! the covars are optimized in the log space
gVarcovsKL = 0.5 - 0.5*vardistx.covars(:)';

gVarmeans = gVarmeans0 + gVarmeans1 + gVarmeans2 + gVarmeansKL;
gVarcovs = gVarcovs0 + gVarcovs1 + gVarcovs2 + gVarcovsKL;

g = [gVarmeans gVarcovs];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
function [gPsi0, gPsi1, gPsi2] = vargpCovGrads(model, vardistx, my, indexPresent)
%
%

d = prod(size(indexPresent));

% change the data (by including the new point and taking only the present indices)
model.m = model.m(:,indexPresent);
model.m = [my; model.m];

pointPsi1 = kernVardistPsi1Compute(model.kern, vardistx, model.X_u);
pointPsi2 = kernVardistPsi2Compute(model.kern, vardistx, model.X_u);

model.Psi1 = [pointPsi1; model.Psi1];

model.C = model.invLm * (model.Psi2 + pointPsi2) * model.invLmT;
model.TrC = sum(diag(model.C)); % Tr(C)
model.At = (1/model.beta) * eye(size(model.C,1)) + model.C;
model.Lat = jitChol(model.At)';
model.invLat = model.Lat\eye(size(model.Lat,1));
model.invLatT = model.invLat';
model.logDetAt = 2*(sum(log(diag(model.Lat)))); % log |At|
model.P1 = model.invLat * model.invLm; % M x M
model.P = model.P1 * (model.Psi1' * model.m);
model.TrPP = sum(sum(model.P .* model.P));
model.B = model.P1' * model.P;
P1TP1 = (model.P1' * model.P1);
Tb = (1/model.beta) * d * P1TP1;
Tb = Tb + (model.B * model.B');
model.T1 = d * model.invK_uu - Tb;

gPsi2 = (model.beta/2) * model.T1;

gPsi0 = -0.5 * model.beta * d;

gPsi1 = model.beta*(P1TP1*model.Psi1'*model.m*my');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
function [gPointDyn, model] = dynPointLogLikeGradient(model, vardistx, y, indexPresent, indexMissing)
jitter = 1e-6;

%%% RE-OPT-CODE-NEW_
if isfield(model.dynamics, 'reoptimise') && model.dynamics.reoptimise
    % In this case the inducing points and theta_t are optimised and
    % since updateStats is not called, any computations which need X_u
    % or theta_t must be done here.
    model = vargplvmDynamicsUpdateStats(model);
    model.K_uu = kernCompute(model.kern, model.X_u);
    model.K_uu = model.K_uu ...
        + sparseDiag(repmat(jitter, size(model.K_uu, 1), 1));
    model.Lm = jitChol(model.K_uu)';
    model.invLm = model.Lm\eye(model.k);
    model.invLmT = model.invLm';
    model.invK_uu = model.invLmT * model.invLm;
end %%% _RE-OPT-CODE-NEW

%%%%%%%%% Only the missing parts and from the training data only
mOrig = model.m;
vardistDynTemp = model.dynamics.vardist;

N = model.N;
Nstar = size(y,1);
Kt = zeros(N+Nstar, N+Nstar);


%model.y = model.y(:,indexMissing);
model.m = model.m(:,indexMissing);
model.d = prod(size(indexMissing));
% Augment the time vector to include the timestamp of the new point
model.dynamics.t = [model.dynamics.t; model.dynamics.t_star];
% Augment the reparametrized variational parameters mubar and lambda
model.dynamics.vardist.means = [model.dynamics.vardist.means; vardistx.means];
model.dynamics.vardist.covars = [model.dynamics.vardist.covars; vardistx.covars];
model.dynamics.N = model.dynamics.N + Nstar;
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
% cross block
if ~isfield(model.dynamics, 'seq') || isempty(model.dynamics.seq)
    Kt(1:N, N+1:end) = kernCompute(model.dynamics.kern, model.dynamics.t(1:N), model.dynamics.t_star);
    Kt(N+1:end, 1:N) = Kt(1:N, N+1:end)';
end
model.dynamics.Kt = Kt;
model.dynamics.fixedKt = 1;
model = vargpTimeDynamicsUpdateStats(model);

vardist2 = model.vardist;
vardist2.means = model.vardist.means(1:end-Nstar,:);
vardist2.covars = model.vardist.covars(1:end-Nstar,:);
vardist2.nParams = 2*prod(size(vardist2.means));
vardist2.numData = size(vardist2.means,1);
model.Psi0 = kernVardistPsi0Compute(model.kern, vardist2);
model.Psi1 = kernVardistPsi1Compute(model.kern, vardist2, model.X_u);
model.Psi2 = kernVardistPsi2Compute(model.kern, vardist2, model.X_u);

model.C = model.invLm * model.Psi2 * model.invLmT;
%model.TrC = sum(diag(model.C)); % Tr(C)
model.At = (1/model.beta) * eye(size(model.C,1)) + model.C; % At = beta^{-1} I + C
model.Lat = jitChol(model.At)';
model.invLat = model.Lat\eye(size(model.Lat,1));
%model.invLatT = model.invLat';
%model.logDetAt = 2*(sum(log(diag(model.Lat)))); % log |At|
model.P1 = model.invLat * model.invLm; % M x M
model.P = model.P1 * (model.Psi1' * model.m);
%model.TrPP = sum(sum(model.P .* model.P));
model.B = model.P1' * model.P;
%model.invK_uu = model.invLmT * model.invLm;
Tb = (1/model.beta) * model.d * (model.P1' * model.P1);
Tb = Tb + (model.B * model.B');
model.T1 = model.d * model.invK_uu - Tb;

modelTmp = model;
modelTmp.vardist = vardist2;
modelTmp.dynamics.vardist = vardistDynTemp;
modelTmp.dynamics.Kt = modelTmp.dynamics.Kt(1:end-Nstar,1:end-Nstar);
modelTmp.dynamics.N = N;
modelTmp.dynamics.t = modelTmp.dynamics.t(1:end-Nstar);
modelTmp.dynamics.vardist.numData = modelTmp.dynamics.N; %%%% new...
modelTmp.vardist.numData = modelTmp.N;%%%% new...
%modelTmp.TrYY = sum(sum(model.m .* model.m));



% GRADIENT OF THE LL1 TERM mu,S from the training data

% To understand the following code it's helpful to see
% vargpTimeDynamicsPriorReparamGrads.m; LL1 term now is not computed
% together with KL so we cannot use that function here; however, it does
% contain theta_t so these computations must be done but only for the part
% corresponding to modelTmp. Thus, the following code actually is a copy of
% vargpTimeDynamicsPriorReparamGrads but is applied only on the appropriate
% part of Kt.
% i.e: gPointDyn1Tmp(:,1:model.dynamics.q) -> gVarmeansLik
%      gPointDyn1Tmp(:,model.dynamics.q+q) -> gVarcovsLik(:,q)

%gPointDyn1 = zeros(Nstar,modelTmp.dynamics.q*2);
gPointDyn1 = zeros(N+Nstar,modelTmp.dynamics.q*2);
if ~isempty(indexMissing)
    if isfield(model.dynamics, 'reoptimise') && model.dynamics.reoptimise %%% RE-OPT-CODE-NEW
        % gDynKern will be empty in that case (includeKL is 0) %%% RE-OPT-CODE-NEW
        [gPointDyn1Tmp, gDynKern, gInd1]= vargplvmLogLikeGradientsVar1(modelTmp, 0); %%% RE-OPT-CODE-NEW
        sumTrGradKL = 0; %%% RE-OPT-CODE-NEW
    else %%% RE-OPT-CODE-NEW
        gPointDyn1Tmp = vargplvmLogLikeGradientsVar1(modelTmp, 0);
    end %%% RE-OPT-CODE-NEW
    gPointDyn1Tmp = reshape(gPointDyn1Tmp, modelTmp.dynamics.vardist.numData, modelTmp.dynamics.q*2);
    
    % means (Note: Kt(1:N,:)' == Kt(:,1:N)  )
    gPointDyn1(:,1:modelTmp.dynamics.q) = model.dynamics.Kt(1:N,:)'*gPointDyn1Tmp(:,1:model.dynamics.q);
    
    % covars
    %gcovTmp = zeros(Nstar,modelTmp.dynamics.q);
    gcovTmp = zeros(N+Nstar,modelTmp.dynamics.q);   
    sumTrGradKL1 = zeros(N+Nstar,N); 
    sumTrGradKL2 = zeros(N+Nstar, N+Nstar); 
    for q=1:model.dynamics.q
        LambdaH_q = model.dynamics.vardist.covars(:,q).^0.5;
        Bt_q = eye(model.dynamics.N) + LambdaH_q*LambdaH_q'.*model.dynamics.Kt;
        % Invert Bt_q
        Lbt_q = jitChol(Bt_q)';
        G1 = Lbt_q \ diag(LambdaH_q);
        G = G1*model.dynamics.Kt;
        % Find Sq
        Sq = model.dynamics.Kt - G'*G;
        Sq = - (Sq .* Sq);
        
        % only the cross matrix is needed as in barmu case
        %gcovTmp(:,q) = Sq(1:N,N+1:end)'*gPointDyn1Tmp(:,model.dynamics.q+q);
        
        % gPointDyn1Tmp contains [means covars]. model.dynamics.q+q gives
        % us only the covars.
        gcovTmp(:,q) = Sq(1:N,:)'*gPointDyn1Tmp(:,model.dynamics.q+q);
        
        %%% RE-OPT-CODE-NEW_   
        if isfield(model.dynamics, 'reoptimise') && model.dynamics.reoptimise
            % Find the coefficient for the grad. wrt theta_t (params of Kt)
            G1T=G1';
            Bhat=G1T*G1;
            BhatKt=G1T*G;
            
                  
            IBK = eye(model.dynamics.N) - BhatKt;

            diagVarcovs = repmat(gPointDyn1Tmp(:,model.dynamics.q+q)', model.dynamics.N,1);         
            sumTrGradKL2 = sumTrGradKL2  +  (IBK(:,1:N).*diagVarcovs)*IBK(:,1:N)';
            
           %%%
         %   tmp = [gPointDyn1Tmp(:,model.dynamics.q+q)' zeros(1,Nstar)]; %tmp
         %   tmp2=IBK*diag(tmp)*IBK'; %tmp
         %   sumTrGradKL2 = sumTrGradKL2  +  tmp2; % tmp
           %%% 

            
            
            sumTrGradKL1 = sumTrGradKL1 + model.dynamics.vardist.means(:,q)*gPointDyn1Tmp(:,q)';
            %if ~isempty(gInd)
            %    trGradKL = trGradKL + dynModel.vardist.means(:,q) * gInd(:,q)';
            %end
        end
        %%% _RE-OPT-CODE-NEW
    end
    
    gPointDyn1(:,(modelTmp.dynamics.q+1):(modelTmp.dynamics.q*2)) = gcovTmp;
    
    if isfield(model.dynamics, 'reoptimise') && model.dynamics.reoptimise                   %%% RE-OPT-CODE-NEW
        gDynKern1 = kernGradient(model.dynamics.kern, model.dynamics.t, model.dynamics.t(1:N),  sumTrGradKL1)  ; %%% RE-OPT-CODE-NEW
        gDynKern1 = gDynKern1 + kernGradient(model.dynamics.kern, model.dynamics.t, sumTrGradKL2); %%% RE-OPT-CODE-NEW
    end                                                                                     %%% RE-OPT-CODE-NEW
    
end


% GRADIENT FOR LL2 plus KL TERMS
%
model.m = mOrig;
% normalize y exactly as model.m is normalized
my = y - repmat(model.bias(indexPresent), Nstar, 1);
my = my./repmat(model.scale(indexPresent), Nstar, 1);
model.m = model.m(:,indexPresent);
model.m = [model.m; my];
model.N = model.N + Nstar;
d = prod(size(indexPresent));
model.d = d;
model.dynamics.nParams = model.dynamics.nParams + 2*prod(size(vardistx.means));
model.nParams = model.nParams + 2*prod(size(vardistx.means));

vardist2 = model.vardist;
vardist2.means = model.vardist.means(end-Nstar+1:end,:);
vardist2.covars = model.vardist.covars(end-Nstar+1:end,:);
vardist2.nParams = 2*prod(size(vardist2.means));
vardist2.numData = size(vardist2.means,1);
pointPsi0 = kernVardistPsi0Compute(model.kern, vardist2);
pointPsi1 = kernVardistPsi1Compute(model.kern, vardist2, model.X_u);
pointPsi2 = kernVardistPsi2Compute(model.kern, vardist2, model.X_u);
model.Psi1 = [model.Psi1; pointPsi1];
model.Psi2 = model.Psi2 + pointPsi2;
model.Psi0 = model.Psi0 + pointPsi0;

model.C = model.invLm * model.Psi2 * model.invLmT;
%model.TrC = sum(diag(model.C)); % Tr(C)
model.At = (1/model.beta) * eye(size(model.C,1)) + model.C; % At = beta^{-1} I + C
model.Lat = jitChol(model.At)';
model.invLat = model.Lat\eye(size(model.Lat,1));
%model.invLatT = model.invLat';
%model.logDetAt = 2*(sum(log(diag(model.Lat)))); % log |At|
model.P1 = model.invLat * model.invLm; % M x M
model.P = model.P1 * (model.Psi1' * model.m);
%model.TrPP = sum(sum(model.P .* model.P));
model.B = model.P1' * model.P;
%model.invK_uu = model.invLmT * model.invLm;
Tb = (1/model.beta) * model.d * (model.P1' * model.P1);
Tb = Tb + (model.B * model.B');
model.T1 = model.d * model.invK_uu - Tb;

if isfield(model.dynamics, 'reoptimise') && model.dynamics.reoptimise %%% RE-OPT-CODE-NEW
    [gPointDyn2, gDynKern2, gInd2] = vargplvmLogLikeGradientsVar1(model,1); %%% RE-OPT-CODE-NEW
else %%% RE-OPT-CODE-NEW
    gPointDyn2 = vargplvmLogLikeGradientsVar1(model,1);
end %%% RE-OPT-CODE-NEW

gPointDyn2 = reshape(gPointDyn2, model.dynamics.vardist.numData, model.dynamics.q*2);
%gPointDyn2 = gPointDyn2(N+1:end,:);

gPointDyn = gPointDyn1 + gPointDyn2;

%%% RE-OPT-CODE-NEW_
gPointDyn = gPointDyn(:)';
if isfield(model.dynamics, 'reoptimise') && model.dynamics.reoptimise
    gDynKern = gDynKern1 + gDynKern2;
    gInd = gInd1 + gInd2;
    if strcmp(model.dynamics.kern.comp{1}.type,'rbf') || strcmp(model.dynamics.kern.comp{1}.type,'matern32')
        if ~isfield(model.dynamics, 'learnVariance') || ~model.dynamics.learnVariance
            gDynKern(2) = 0;
        end
    end
    gPointDyn = [gPointDyn 0*gDynKern gInd]; %%%%%%%%%% TEMP!!!!
end
%%% _RE-OPT-CODE-NEW



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
%function gVar = vargplvmLogLikeGradientsVar1(model, includeKL) %%% RE-OPT-CODE-REM
function [gVar, gDynKern, gInd] = vargplvmLogLikeGradientsVar1(model, includeKL) %%% RE-OPT-CODE-NEW
% Like vargplvmLogLikeGradients but only for the variational parameters (ehm, not exactly). %%% RE-OPT-CODE-MOD


% Likelihood terms (coefficients)
gPsi1 = model.beta * model.m * model.B';
gPsi1 = gPsi1'; % because it is passed to "kernVardistPsi1Gradient" as gPsi1'...
gPsi2 = (model.beta/2) * model.T1;
gPsi0 = -0.5 * model.beta * model.d;

[gKern1, gVarmeans1, gVarcovs1, gInd1] = kernVardistPsi1Gradient(model.kern, model.vardist, model.X_u, gPsi1');
[gKern2, gVarmeans2, gVarcovs2, gInd2] = kernVardistPsi2Gradient(model.kern, model.vardist, model.X_u, gPsi2);
[gKern0, gVarmeans0, gVarcovs0] = kernVardistPsi0Gradient(model.kern, model.vardist, gPsi0);

gVarmeansLik = gVarmeans0 + gVarmeans1 + gVarmeans2;

gVarcovsLik = gVarcovs0 + gVarcovs1 + gVarcovs2;
if includeKL
    %[gVarmeans gVarcovs gDynKern] = modelPriorReparamGrads(model.dynamics, gVarmeansLik, gVarcovsLik);
    % TODO: This may have to be changed if we have model.fixInducing==1.
    [gVarmeans gVarcovs gDynKern] = modelPriorReparamGrads(model.dynamics, gVarmeansLik, gVarcovsLik,[]);
else
    dynModel = model.dynamics;
    gVarmeansLik = reshape(gVarmeansLik,dynModel.N,dynModel.q);
    gVarcovsLik = reshape(gVarcovsLik,dynModel.N,dynModel.q);
    gVarmeans = gVarmeansLik;
    gVarcovs = gVarcovsLik;
    
    gVarcovs = gVarcovs(:)';
    gVarmeans = gVarmeans(:)';
end

gVarcovs = (gVarcovs(:).*model.dynamics.vardist.covars(:))';
gVar = [gVarmeans gVarcovs];

%%% RE-OPT-CODE-NEW_
if isfield(model.dynamics, 'reoptimise') && model.dynamics.reoptimise
    % Compute Gradients with respect to X_u
    gK_uu = 0.5 * (model.T1 - (model.beta * model.d) * model.invLmT * model.C * model.invLm);
    gKX = kernGradX(model.kern, model.X_u, model.X_u);
    gKX = gKX*2;
    dgKX = kernDiagGradX(model.kern, model.X_u);
    for i = 1:model.k
        gKX(i, :, i) = dgKX(i, :);
    end
    gX_u = zeros(model.k, model.q);
    for i = 1:model.k
        for j = 1:model.q
            gX_u(i, j) = gKX(:, j, i)'*gK_uu(:, i);
        end
    end
    gInd = gInd1 + gInd2 + gX_u(:)';
    if ~includeKL
        gDynKern = [];
    end
end
%%% _RE-OPT-CODE-NEW

