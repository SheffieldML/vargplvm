function ll = vargplvmLogLikelihood(model)

% VARGPLVMLOGLIKELIHOOD Log-likelihood for a variational GP-LVM.
% FORMAT
% DESC returns the log likelihood for a given GP-LVM model.
% ARG model : the model for which the log likelihood is to be
% computed. The model contains the data for which the likelihood is
% being computed in the 'y' component of the structure.
% RETURN ll : the log likelihood of the data given the model.
%
% COPYRIGHT : Michalis K. Titsias, 2009-2011
% COPYRIGHT : Neil D. Lawrence, 2009-2011
% COPYRIGHT : Andreas Damianou, 2010-2011


% VARGPLVM

% Note: The 'onlyKL' and 'onlyLikelihood' fields can be set by external
% wrappers and cause the function to only calculate the likelihood or the
% KL part of the variational bound.

% Likelihood term
if ~(isfield(model, 'onlyKL') && model.onlyKL)
    if length(model.beta)==1
        ll = -0.5*(model.d*(-(model.N-model.k)*log(model.beta) ...
            + model.logDetAt) ...
            - (model.TrPP ...
            - model.TrYY)*model.beta);
        %if strcmp(model.approx, 'dtcvar')
        ll = ll - 0.5*model.beta*model.d*model.Psi0 + 0.5*model.d*model.beta*model.TrC;
        %end
    else
        error('Not implemented variable length beta yet.');
    end
    ll = ll-model.d*model.N/2*log(2*pi);
else
    ll=0;
end

%{
% %%%%%%%%TEMP_ (for analysing the "difficult" trace term)
%try
%    load TEMPLikelihoodTrace
%catch exception
%    TEMPLikelihoodTrace=[];
%end
%TEMPLikelihoodTrace=[TEMPLikelihoodTrace ;-0.5*model.beta*model.d*model.Psi0 + 0.5*model.d*model.beta*model.TrC];
%save 'TEMPLikelihoodTrace.mat' 'TEMPLikelihoodTrace';
%%% _TEMP
% load TEMPbetaLikTrC;
% TEMPbetaLikTrC = [TEMPbetaLikTrC 0.5*model.d*model.beta*model.TrC];
% save 'TEMPbetaLikTrC.mat' 'TEMPbetaLikTrC';
%
% load TEMPbetaLikNksigm;
% TEMPbetaLikNksigm=[TEMPbetaLikNksigm (model.d*(-(model.N-model.k)*log(model.beta)))];
% save 'TEMPbetaLikNksigm.mat' 'TEMPbetaLikNksigm';
%
% load TEMPbetaLikPsi0;
% TEMPbetaLikPsi0=[TEMPbetaLikPsi0 (- 0.5*model.beta*model.d*model.Psi0)];
% save 'TEMPbetaLikPsi0.mat' 'TEMPbetaLikPsi0';
%
% load TEMPbetaLikTrPP;
% TEMPbetaLikTrPP=[TEMPbetaLikTrPP 0.5*model.TrPP];
% save 'TEMPbetaLikTrPP.mat' 'TEMPbetaLikTrPP';
%
% load TEMPbetaLikLat;
% TEMPbetaLikLat=[TEMPbetaLikLat -0.5*model.d*model.logDetAt];
% save 'TEMPbetaLikLat.mat' 'TEMPbetaLikLat';
%%%%%%%%%%%%
%}

% KL divergence term
if ~(isfield(model, 'onlyLikelihood') && model.onlyLikelihood)
    if isfield(model, 'dynamics') && ~isempty(model.dynamics)
        % A dynamics model is being used. (The following is actually
        % the -KL term)
        KLdiv = modelVarPriorBound(model);
        KLdiv = KLdiv + 0.5*model.q*model.N; %%% The constant term!!
    else
        varmeans = sum(sum(model.vardist.means.*model.vardist.means));
        varcovs = sum(sum(model.vardist.covars - log(model.vardist.covars)));
        KLdiv = -0.5*(varmeans + varcovs) + 0.5*model.q*model.N;
    end
else
    KLdiv=0;
end

% Obtain the final value of the bound by adding the likelihood
% and the KL term.
ll = ll + KLdiv;
