function model = vargplvmAddParamPrior(model, paramName, priorName, varargin)
% VARGPLVMADDPARAMPRIOR Add a prior in one of the parameters of the given
% vargplvm model
%
% DESC Add a prior in one of the parameters of the model.  This prior will
%      NOT be optimised (ie fixed parameters).
%
% ARG model: The model for which we want to add a prior
% ARG paramName: The name of the parameter(s) for which the prior will be
%                added (a regexp will match the given name)
% ARG priorName: The prior to be used. Make sure that the function
%                [priorName 'priorParamInit'] exists.
% ARG varargin: Any arguments needed for the specific constructor of the
%               prior are going to be pased
% 
% EXAMPLE: 
% % Add a prior on the noise parameter to encourage high SNR. High SNR is
% % achieved when beta is high compared to the data variance. So we will
% % constrain it relative to that value.
% >> varData = var(model.m(:));
% >> meanP = 400; stdP = 100; % Plot this with: Nsamp = 10000;XX = zeros(Nsamp,1);for iii=1:Nsamp, XX(iii,:) = meanP + randn.*stdP;end; hist(XX,100);
% % Mean and std of the prior depend on the variance of the data, this says
% % that no matter what the variance of the data is, the SNR is encouraged to
% % be around meanP. 
% >> prec = 1/(stdP^2*varData);
% >> mu = meanP/varData;
% % prec is actually multiplied by the gradient. Should be selected so that
% % the overall gradient coming from this precision is comparable to the
% % gradient without the prior.
% >> model = vargplvmAddParamPrior(model, 'beta', 'gaussian2', [prec mu]);
% % Optional, for stronger effect. For beta, a reasonable choice is model.N, because
% % beta is actually N parameters tied together (e.g. if we had
% % heteroscedastic noise) and it appears as a sum over N everywhere in the
% % likelihood term, so we have to do the same in the prior p(beta).
% % !!!!!!!!!!! TODO: This might be wrong, because we compute the log
% % probability
% >> model.paramPriors{1}.prior.scale = model.N;
%
% % For example, for the gaussian2 prior, the gradient is:
% g = -(prior.precision*length(x)).*(x - prior.mean);
% so we select precision to be:
% g = vargplvmLogLikeGradients(model);
% ch = transformedParamChainRule(model.betaTransform, 1, model.beta);
% prec = g(model.paramPriors{1}.index)/(ch*model.paramPriors{1}.prior.dim*(model.beta-mu))
% model.paramPriors{1}.prior.precisions = abs(prec);
%
% EXAMPLE 2: Inverse Gamma prior - TODO (there are mistakes)!!!!!!!!!!!
%{
%varData = var(model.m(:));
meanP = 100; % The peak of the distribution for the SNR
Nsamp = 10000;
XX = zeros(Nsamp,1);
for j=1:Nsamp
    XX(j,:) = 1./varData * invgamrnd(4.2,meanP*5.7);
end
subplot(1,2,1);
hist(XX,500); title('distribution for beta prior')
subplot(1,2,2);
hist(XX.*varData,500); title('distribution for SNR prior')

%%

SNR_range = 0:0.5:2000; beta_range = SNR_range;
prior.a = 4.2; prior.b = meanP*5.7;
ll = zeros(length(beta_range),1);
for j=1:length(beta_range)
    ll(j,:) = log(1/varData) + invgammaPriorLogProb(prior, beta_range(j)); % Distribution for beta
end
subplot(1,2,1);
plot(beta_range,exp(ll)) ;title('distribution for beta prior')
subplot(1,2,2);
plot(SNR_range,exp(ll.*varData)) ; title('distribution for SNR prior')

%}
%
% SEEALSO: vargplvmParamPriorLogProb.m, vargplvmParamPriorGradients.m,
% priorParamInit.m
% 
% COPYRIGHT: Andreas C. Damianou, 2013
%
% VARGPLVM

[~, names] = vargplvmExtractParam(model);

if ~isfield(model, 'paramPriors')
    model.paramPriors = {};
end

%--- INDEX
onlyFirst = false;
switch paramName
    case {'beta', 'input scale'}
        tmpIndex = cellfun(@(x)regexpi(x,paramName),names,'un',0);
    case {'dynamicsWhiteKernelVariance'}
        searchName = 'Kernel, white 1 variance';
        tmpIndex = cellfun(@(x)regexpi(x,searchName),names,'un',0);
        onlyFirst = true;
    otherwise
        error('Prior requested for unrecognised parameter')
end

paramPriors.index = [];
for i=1:length(tmpIndex)
    if ~isempty(tmpIndex{i})
        paramPriors.index = [paramPriors.index i];
        if onlyFirst
            break
        end
    end
end
fprintf('# Added a %s prior on parameter: %s (%d indices)\n', priorName, paramName, length(paramPriors.index));

%---- PARAM. NAME
paramPriors.paramName = paramName;

%--- Actual prior
paramPriors.prior = priorCreate(priorName, varargin{:});

% Add prior
model.paramPriors{end+1} = paramPriors;

