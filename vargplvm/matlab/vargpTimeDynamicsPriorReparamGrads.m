function [gVarmeans gVarcovs gDynKern] = vargpTimeDynamicsPriorReparamGrads(dynModel, gVarmeansLik, gVarcovsLik, gInd)
% VARGPTIMEDYNAMICSPRIORREPARAMGRADS Returns the gradients of the various types of the
% variational GPLVM bound when dynamics is used. The gradients returned are
% those for the variational means and covariances and the temporal kernel.
% FORMAT
% DESC The function takes a dynamics model structure along with some partial
% derivatives, amends as appropriate  and calculates the whole gradient for the
% variational bound. The gradients returned are w.r.t mu_bar and lambda, i.e. the parameters used by the optimiser
% (the ones introduced by the reparametrization), not the original ones (mu and S).
% These gradients are obtained in two steps: firstly, another funtion is used to calculate the
% gradients that correspond only to the likelihood part of the bound and
% are w.r.t the original parameters mu and S. These quantities are the arguments
% gVarmenasLik and gVarcovsLik and are calculated using the already
% implemented code for the static vargplvm, which assumes that the original
% parameters are not coupled.
% The current funtions receives these quantities and a) amends
% with partial derivatives because the variational parameters are coupled
% via Kt for the dyn. gplvm b) Computes the whole gradients for both parts
% of the bound, the one corresponding to the likelihood and the one
% corresponding to the prior. This must be done in a single function
% because the final formula contains both parts in a nonlinear form.
% c) Also the kernel hyperparameters for the dynamics kernel are being retured.
% See the dyn. vargplvm notes for more details.
%
% ARG dynModel : the dynamics model structure for which the gradients are
% to be computed.
% ARG gVarmeansLik, gVarcovsLik: the gradients for the VAR-GPLVM model computed
% for the original parameters and only for the likelihood term.
% ARG gInd: in case the inducing points are tied to the variational means
% this is the partial derivatives of the likelihood part of the variational
% bound w.r.t the inducing points X_u, otherwise it is just [].
% RETURN gVarmeans, gVarcovs : the gradients for the "reparametrized" means
% and covariances (the ones that are visible to the optimiser, not the
% original ones) for the VAR-GPLVM model.
% RETURN gDynKern: the gradient w.r.t the hyperparameters of the dynamics
% kernel
%
% SEEALSO : vargplvmLogLikeGradients modelVarPriorBound
%
% COPYRIGHT : Michalis K. Titsias and Andreas C. Damianou, 2010-2011


% VARGPLVM


if isfield(dynModel, 'seq') & ~isempty(dynModel.seq)
    [gVarmeans gVarcovs gDynKern] = vargpTimeDynamicsPriorReparamGradsSeq(dynModel, gVarmeansLik, gVarcovsLik, gInd);
else
    [gVarmeans gVarcovs gDynKern] = vargpTimeDynamicsPriorReparamGrads1(dynModel, gVarmeansLik, gVarcovsLik, gInd);
end


%-------------
% The following function calculates the derivative for the bound term
% corresponding to the prior, when the model learns from multiple
% independent sequences.
function [gVarmeans gVarcovs gDynKern] = vargpTimeDynamicsPriorReparamGradsSeq(dynModel, gVarmeansLik, gVarcovsLik, gInd)
% gVarmeansLik and gVarcovsLik are serialized into 1x(NxQ) vectors.
% Convert them to NxQ matrices, i.e. fold them column-wise.
gVarmeansLik = reshape(gVarmeansLik,dynModel.N,dynModel.q);
gVarcovsLik = reshape(gVarcovsLik,dynModel.N,dynModel.q);

% In case gInd is empty then the inducing points are not reparametrized
% (are not fixed to the variational means).
if ~isempty(gInd)
    gInd = reshape(gInd, dynModel.N, dynModel.q); %%%%%
end

gVarmeans = dynModel.Kt * (gVarmeansLik - dynModel.vardist.means);

gVarcovs = zeros(dynModel.N, dynModel.q); % memory preallocation
seqStart=1;
seq = dynModel.seq;
for i=1:length(dynModel.seq)
    seqEnd = seq(i);
    sumTrGradKL{i} = zeros(seqEnd-seqStart+1, seqEnd-seqStart+1);
    seqStart = seqEnd+1;
end


for q=1:dynModel.q
    LambdaH_q = dynModel.vardist.covars(:,q).^0.5;
    
    % The calculations are performed for each sequence independently
    % because correlations between sequences are not captured and, thus,
    % most of the matrices involved are block diagonal, so that in each of
    % the following loops only one block is considered.
    seqStart=1;
    for i=1:length(dynModel.seq)
        seqEnd = seq(i);
        Bt_q = eye(seqEnd-seqStart+1) + LambdaH_q(seqStart:seqEnd,1)*LambdaH_q(seqStart:seqEnd,1)'.*dynModel.Kt(seqStart:seqEnd,seqStart:seqEnd);
        
        % Invert Bt_q
        Lbt_q = jitChol(Bt_q)';
        G1 = Lbt_q \ diag(LambdaH_q(seqStart:seqEnd,1));
        G = G1*dynModel.Kt(seqStart:seqEnd, seqStart:seqEnd);
        
        % Find Sq
        Sq = dynModel.Kt(seqStart:seqEnd, seqStart:seqEnd) - G'*G;
        
        gVarcovs(seqStart:seqEnd,q) = - (Sq .* Sq) * (gVarcovsLik(seqStart:seqEnd,q) + 0.5*dynModel.vardist.covars(seqStart:seqEnd,q));
        
        % Find the coefficient for the grad. wrt theta_t (params of Kt)
        G1T=G1';
        Bhat=G1T*G1;
        BhatKt=G1T*G;
        
        trGradKL = -0.5*(BhatKt*Bhat + dynModel.vardist.means(seqStart:seqEnd,q) * dynModel.vardist.means(seqStart:seqEnd,q)');
        
        IBK = eye(seqEnd-seqStart+1) - BhatKt;
        diagVarcovs = repmat(gVarcovsLik(seqStart:seqEnd,q)', seqEnd-seqStart+1,1);
        trGradKL = trGradKL + IBK .*  diagVarcovs * IBK';
        trGradKL = trGradKL + dynModel.vardist.means(seqStart:seqEnd,q) * gVarmeansLik(seqStart:seqEnd,q)';
        
        % In case gInd is empty then the inducing points are not reparametrized
        % (are not fixed to the variational means) we need not amend further the
        % derivative w.r.t theta_t.
        if ~isempty(gInd)
            trGradKL = trGradKL + dynModel.vardist.means(seqStart:seqEnd,q) * gInd(seqStart:seqEnd,q)';
        end
              
                
        %tmp = eye(dynModel.N);
        %tmp(seqStart:seqEnd, seqStart:seqEnd) = trGradKL;
        %        sumTrGradKL = sumTrGradKL + trGradKL;
        sumTrGradKL{i} = sumTrGradKL{i} + trGradKL;
        % size(sumTrGradKL{i} )
        %sumTrGradKL = sumTrGradKL + tmp;
        
        seqStart = seqEnd+1;
    end
    %sumTrGradKL = sumTrGradKL + trGradKL;
end

seqStart=1;
gDynKern = size(dynModel.kern.nParams,1);
for i=1:length(dynModel.seq)
    seqEnd = seq(i);
    gDynKern = gDynKern + kernGradient(dynModel.kern, dynModel.t(seqStart:seqEnd), sumTrGradKL{i});
    seqStart = seqEnd+1;
end

% Serialize (unfold column-wise) gVarmeans and gVarcovs from NxQ matrices
% to 1x(NxQ) vectors
gVarcovs = gVarcovs(:)';
gVarmeans = gVarmeans(:)';



function [gVarmeans gVarcovs gDynKern] = vargpTimeDynamicsPriorReparamGrads1(dynModel, gVarmeansLik, gVarcovsLik, gInd)

% gVarmeansLik and gVarcovsLik are serialized into 1x(NxQ) vectors.
% Convert them to NxQ matrices, i.e. fold them column-wise.
gVarmeansLik = reshape(gVarmeansLik,dynModel.N,dynModel.q);
gVarcovsLik = reshape(gVarcovsLik,dynModel.N,dynModel.q);

if ~isempty(gInd)
    gInd = reshape(gInd, dynModel.N, dynModel.q); %%%%%
end

% Flag, indicating/multiplier that we need derivatives only w.r.t the
% likelihood part of the bound.
if ~(isfield(dynModel, 'onlyLikelihood') && dynModel.onlyLikelihood)
    onlyLikelihood = 0;
else
    onlyLikelihood = 1;
end

% The second term in the parenthesis, corresponds to the KL term. The
% multiplier will set it to zero, if we need only the derivatives
% corresponding only to the likelihood part.
gVarmeans = dynModel.Kt * (gVarmeansLik - (~onlyLikelihood) * dynModel.vardist.means);

gVarcovs = zeros(dynModel.N, dynModel.q); % memory preallocation
sumTrGradKL = 0;

for q=1:dynModel.q
    LambdaH_q = dynModel.vardist.covars(:,q).^0.5;
    Bt_q = eye(dynModel.N) + LambdaH_q*LambdaH_q'.*dynModel.Kt;
    
    % Invert Bt_q
    Lbt_q = jitChol(Bt_q)';
    G1 = Lbt_q \ diag(LambdaH_q);
    G = G1*dynModel.Kt;
    
    % Find Sq
    Sq = dynModel.Kt - G'*G;
    
    % If onlyLikelihood is set to 1, then the KL part will be multiplied
    % with zero, thus being ignored.
    gVarcovs(:,q) = - (Sq .* Sq) * (gVarcovsLik(:,q) + ...
        (~onlyLikelihood)*0.5*dynModel.vardist.covars(:,q));
    
    % Find the coefficient for the grad. wrt theta_t (params of Kt)
    G1T=G1';
    Bhat=G1T*G1;
    BhatKt=G1T*G;
    
    % If we only need the derivatives from the likelihood part of the
    % bound, set the corrsponding KL part to zero.
    if onlyLikelihood
        trGradKL = 0;
    else
        trGradKL = -0.5*(BhatKt*Bhat + dynModel.vardist.means(:,q) * dynModel.vardist.means(:,q)');
    end
    
    IBK = eye(dynModel.N) - BhatKt;
    diagVarcovs = repmat(gVarcovsLik(:,q)', dynModel.N,1);
    trGradKL = trGradKL + IBK .*  diagVarcovs * IBK';
    trGradKL = trGradKL + dynModel.vardist.means(:,q) * gVarmeansLik(:,q)';
    
    % In case gInd is empty then the inducing points are not reparametrized
    % (are not fixed to the variational means) we need not amend further the
    % derivative w.r.t theta_t, otherwise we have to do that.
    if ~isempty(gInd)
        trGradKL = trGradKL + dynModel.vardist.means(:,q) * gInd(:,q)';
    end
    
    
    sumTrGradKL = sumTrGradKL + trGradKL;
end
gDynKern = kernGradient(dynModel.kern, dynModel.t, sumTrGradKL);

% Serialize (unfold column-wise) gVarmeans and gVarcovs from NxQ matrices
% to 1x(NxQ) vectors
gVarcovs = gVarcovs(:)';
gVarmeans = gVarmeans(:)';



%{
%%% Alternative implementations
function [gVarmeans gVarcovs gDynKern] = vargpTimeDynamicsPriorReparamGrads2(dynModel, gVarmeansLik, gVarcovsLik)

%fprintf(1,'vargpTimeDynamicsVarReparamGrads1..\n');%%% DEBUG

% gVarmeansLik and gVarcovsLik are serialized into 1x(NxQ) vectors.
% Convert them to NxQ matrices, i.e. fold them column-wise.
gVarmeansLik = reshape(gVarmeansLik,dynModel.N,dynModel.q);
gVarcovsLik = reshape(gVarcovsLik,dynModel.N,dynModel.q);

%%% TODO: check that the following is true for each q (compare with report)
gVarmeans = dynModel.Kt * (gVarmeansLik - dynModel.vardist.means);

gVarcovs = zeros(dynModel.N, dynModel.q); % memory preallocation
% gDynKern = 0;
 gMeanPart = 0;
sumTrGradKL = 0;


for q=1:dynModel.q

    
    LambdaH_q = dynModel.vardist.covars(:,q).^0.5;
    Bt_q = eye(dynModel.N) + LambdaH_q*LambdaH_q'.*dynModel.Kt;
    
    % Invert Bt_q
    Lbt_q = jitChol(Bt_q)';
    G1 = Lbt_q \ diag(LambdaH_q);
    G = G1*dynModel.Kt;
   
    % Find Sq
    Sq = dynModel.Kt - G'*G;
    
    gVarcovs(:,q) = - (Sq .* Sq) * (gVarcovsLik(:,q) + 0.5*dynModel.vardist.covars(:,q));
  
    
    Bhat = diag(LambdaH_q) * inv(Bt_q) * diag(LambdaH_q);
    BhatKt = Bhat * dynModel.Kt;
    trGradKL = -0.5*(BhatKt*Bhat + dynModel.vardist.means(:,q) * dynModel.vardist.means(:,q)');
   
    %%%%%%% 1st way %%%
    %IBK = eye(dynModel.N) - BhatKt;
    %trGradKL = trGradKL + IBK' * diag(gVarcovsLik(:,q)) * IBK;
    %------------------
    
    %%%%%%% 2nd way %%%%
    Lt=jitChol(dynModel.Kt)'; invLt = Lt \ eye(dynModel.N); invKt = invLt' * invLt;
    %trGradKL = trGradKL + (Sq * invKt) * diag(gVarcovsLik(:,q)) * (invKt*Sq);
    trGradKL = trGradKL + (invKt*Sq) * diag(gVarcovsLik(:,q)) * (Sq*invKt);
    %trGradKL = trGradKL + (Sq * invKt) * diag(gVarcovsLik(:,q)) * (invKt*Sq);

    %------------------
   
    trGradKL = trGradKL + dynModel.vardist.means(:,q) * gVarmeansLik(:,q)';
     sumTrGradKL = sumTrGradKL + trGradKL;
%     %gMeanPart = gMeanPart + dynModel.vardist.means(:,q) * gVarmeansLik(:,q)';
end
 gDynKern = kernGradient(dynModel.kern, dynModel.t, sumTrGradKL);
 %gKernMu = dynModel.vardist.means * gVarmeansLik';%%%
 %gDynKern = gDynKern + kernGradient(dynModel.kern, dynModel.t, gMeanPart); %%%
 


    %gKern = kernGradient(dynModel.kern, dynModel.t, trGradKL);
    %gDynKern = gDynKern + gKern;
%%%%%%%%%%%%%%%%%
    
%     G1T=G1';
%     A1=G1T*G1;
%     A2=G1T*G;
%
%     %-- Calculate the trace term-by-term
%     trGradKL = -0.5*(A2*A1 + dynModel.vardist.means(:,q) * dynModel.vardist.means(:,q)');
%     IBK = eye(dynModel.N) - A2;
%     diagVarcovs = repmat(gVarcovsLik(:,q)', dynModel.N,1);
%     trGradKL = trGradKL + IBK' .*  diagVarcovs * IBK;
%    %     trGradKL = trGradKL + IBK' * diag(gVarcovsLik(:,q)) * IBK;
%     gKern = kernGradient(dynModel.kern, dynModel.t, trGradKL);
%     gDynKern = gDynKern + gKern;
%
%     %-- Calculate the means term
%     % The following might be the other way around (because
%     % kernGradient traces the result -it should be correct even in that case thought-)
%     gMeanPart = gMeanPart + gVarmeansLik(:,q)' * dynModel.vardist.means(:,q);
% end


%gDynKern = gDynKern + kernGradient(dynModel.kern, dynModel.t, gMeanPart);

% The following should rather go to vargplvmLogLikeGradients
%gVarcovs = gVarcovs.*dynModel.vardist.covars;

% Serialize (unfold column-wise) gVarmeans and gVarcovs from NxQ matrices
% to 1x(NxQ) vectors
gVarcovs = gVarcovs(:)';
gVarmeans = gVarmeans(:)';

% From the old one with a small addition
function [gVarmeans gVarcovs gDynKern] = vargpTimeDynamicsPriorReparamGradsOLD(dynModel, gVarmeansLik, gVarcovsLik)


%fprintf(1,'vargpTimeDynamicsVarReparamGrads2..\n');%%% DEBUG

% gVarmeansLik and gVarcovsLik are serialized into 1x(NxQ) vectors.
% Convert them to NxQ matrices, i.e. fold them column-wise.
gVarmeansLik = reshape(gVarmeansLik,dynModel.N,dynModel.q);
gVarcovsLik = reshape(gVarcovsLik,dynModel.N,dynModel.q);

%%% TODO: check that the following is true for each q (compare with report)
gVarmeans = dynModel.Kt * (gVarmeansLik - dynModel.vardist.means);

gVarcovs = zeros(dynModel.N, dynModel.q); % memory preallocation
gDynKern = 0;
gMeanPart = 0;
for q=1:dynModel.q
%    LambdaH_q = diag(dynModel.vardist.covars(:,q))^(0.5); % SLOW!
    LambdaH_q = dynModel.vardist.covars(:,q).^0.5;
    %Bt_q = eye(dynModel.N) + LambdaH_q * dynModel.Kt * LambdaH_q;
    Bt_q = eye(dynModel.N) + LambdaH_q*LambdaH_q'.*dynModel.Kt;
    % Invert Bt_q
    Lbt_q = jitChol(Bt_q)';
    LbtInv_q = Lbt_q \ eye(dynModel.N);
    BtInv_q = LbtInv_q' * LbtInv_q;
    
    % Bhat_q = LambdaH_q * BtInv_q * LambdaH_q;   %SLOW
    
    % Since LambdaH_q is diagonal, we can use the property
    % D*M*D=diag(D)*diag(D)'.*M, D diag.matrix, M symmetric.
    Bhat_q = LambdaH_q*LambdaH_q'.*BtInv_q;
    
    
    BhatKt = Bhat_q * dynModel.Kt;
    % Find Sq
    Sq = dynModel.Kt - dynModel.Kt * BhatKt;
  % gVarcovs(:,q) = (Sq .* Sq) * (gVarcovsLik(:,q) - 0.5*dynModel.vardist.covars(:,q));
    gVarcovs(:,q) = - (Sq .* Sq) * (gVarcovsLik(:,q) + 0.5*dynModel.vardist.covars(:,q));
    
    
    
    %-- Calculate the trace term-by-term
    trGradKL = -0.5*(BhatKt * Bhat_q + dynModel.vardist.means(:,q) * dynModel.vardist.means(:,q)');
    IBK = eye(dynModel.N) - BhatKt;
    
    %%% NEW
    diagVarcovs = repmat(gVarcovsLik(:,q)', dynModel.N,1);
    trGradKL = trGradKL + IBK' .*  diagVarcovs * IBK;
    %%%
    
    %trGradKL = trGradKL + IBK' * diag(gVarcovsLik(:,q)) * IBK;
    gKern = kernGradient(dynModel.kern, dynModel.t, trGradKL);
    gDynKern = gDynKern + gKern;
    
    %-- Calculate the means term
    % The following might be the other way around (because
    % kernGradient traces the result -it should be correct even in that case thought-)
    gMeanPart = gMeanPart + gVarmeansLik(:,q)' * dynModel.vardist.means(:,q);
    %gMeanPart = gVarmeansLik(:,q)' * dynModel.vardist.means(:,q);
    %gDynKern = gDynKern + kernGradient(dynModel.kern, dynModel.t, gMeanPart);
end
gDynKern = gDynKern + kernGradient(dynModel.kern, dynModel.t, gMeanPart);


% Serialize (unfold column-wise) gVarmeans and gVarcovs from NxQ matrices to 1x(NxQ) vectors
gVarcovs = gVarcovs(:)';
gVarmeans = gVarmeans(:)';

%-------------
function [gVarmeans gVarcovs gDynKern] = vargpTimeDynamicsPriorReparamGrads1OLD(dynModel, gVarmeansLik, gVarcovsLik)

% gVarmeansLik and gVarcovsLik are serialized into 1x(NxQ) vectors.
% Convert them to NxQ matrices, i.e. fold them column-wise.
gVarmeansLik = reshape(gVarmeansLik,dynModel.N,dynModel.q);
gVarcovsLik = reshape(gVarcovsLik,dynModel.N,dynModel.q);

gVarmeans = dynModel.Kt * (gVarmeansLik - dynModel.vardist.means);

gVarcovs = zeros(dynModel.N, dynModel.q); % memory preallocation
sumTrGradKL = 0;

for q=1:dynModel.q
    LambdaH_q = dynModel.vardist.covars(:,q).^0.5;
    Bt_q = eye(dynModel.N) + LambdaH_q*LambdaH_q'.*dynModel.Kt;
    
    % Invert Bt_q
    Lbt_q = jitChol(Bt_q)';
    G1 = Lbt_q \ diag(LambdaH_q);
    G = G1*dynModel.Kt;
   
    % Find Sq
    Sq = dynModel.Kt - G'*G;
    
    gVarcovs(:,q) = - (Sq .* Sq) * (gVarcovsLik(:,q) + 0.5*dynModel.vardist.covars(:,q));

    G1T=G1';
    Bhat=G1T*G1;
    BhatKt=G1T*G;
    trGradKL = -0.5*(BhatKt*Bhat + dynModel.vardist.means(:,q) * dynModel.vardist.means(:,q)');
    IBK = eye(dynModel.N) - BhatKt;
   % trGradKL = trGradKL + IBK * diag(gVarcovsLik(:,q)) * IBK'; % Correct
    diagVarcovs = repmat(gVarcovsLik(:,q)', dynModel.N,1);
    trGradKL = trGradKL + IBK .*  diagVarcovs * IBK';
%    Bhat = diag(LambdaH_q) * inv(Bt_q) * diag(LambdaH_q);
%    BhatKt = Bhat * dynModel.Kt;
%    trGradKL = -0.5*(BhatKt*Bhat + dynModel.vardist.means(:,q) * dynModel.vardist.means(:,q)');
%
    %%%%%%% 1st way %%%
    %IBK = eye(dynModel.N) - BhatKt;
    %trGradKL = trGradKL + IBK' * diag(gVarcovsLik(:,q)) * IBK;
    %------------------
    
%     %%%%%%% 2nd way %%%%
%     Lt=jitChol(dynModel.Kt)'; invLt = Lt \ eye(dynModel.N); invKt = invLt' * invLt;
%     %trGradKL = trGradKL + (Sq * invKt) * diag(gVarcovsLik(:,q)) * (invKt*Sq);
%     %trGradKL = trGradKL + (invKt*Sq) * diag(gVarcovsLik(:,q)) * (Sq*invKt); % Correct
%     invKtSq = invKt*Sq;
%     trGradKL = trGradKL + invKtSq * diag(gVarcovsLik(:,q)) * invKtSq'; % Correct
%     %trGradKL = trGradKL + (Sq * invKt) * diag(gVarcovsLik(:,q)) * (invKt*Sq);
%
%     %------------------
   
     trGradKL = trGradKL + dynModel.vardist.means(:,q) * gVarmeansLik(:,q)';
     sumTrGradKL = sumTrGradKL + trGradKL;
end
 gDynKern = kernGradient(dynModel.kern, dynModel.t, sumTrGradKL);

% The following should rather go to vargplvmLogLikeGradients
%gVarcovs = gVarcovs.*dynModel.vardist.covars;

% Serialize (unfold column-wise) gVarmeans and gVarcovs from NxQ matrices
% to 1x(NxQ) vectors
gVarcovs = gVarcovs(:)';
gVarmeans = gVarmeans(:)';
%}
