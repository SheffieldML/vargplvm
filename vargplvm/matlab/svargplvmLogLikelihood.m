function ll = svargplvmLogLikelihood(model)

% SVARGPLVMLOGLIKELIHOOD Log-likelihood for a shared variational GP-LVM.
% FORMAT
% DESC returns the log likelihood for a given SVARGP-LVM model.
% ARG model : the model for which the log likelihood is to be
% computed. The model contains the data for which the likelihood is
% being computed in the 'y' component of the structure.
% RETURN ll : the log likelihood of the data given the model.
%
% SEEALSO:  svargplvmLogLikeGradients
%
% COPYRIGHT : Andreas C. Damianou, 2011

% VARGPLVM

% Functions f1 and f2 should be equivalent in the static case, but f1 might
% (?) be faster.
if ~isfield(model, 'dynamics') || isempty(model.dynamics)
    ll = f1(model);
    %ll = f2(model);
else
    ll = f2(model);
    for i = 1:model.numModels
        if isfield(model.comp{i}, 'KLweight')
            assert(model.comp{i}.KLweight == 0.5); % not implemented yet for dynamics
        end
    end
end



function f = f1(model)

% model.numModels=1; %%%%%%%%TEMP

% This works only when there are NOT any dynamics:
varmeans = sum(sum(model.vardist.means.*model.vardist.means));
varcovs = sum(sum(model.vardist.covars - log(model.vardist.covars)));
KLdiv = -0.5*(varmeans + varcovs) + 0.5*model.q*model.N;

model = svargplvmPropagateField(model, 'onlyLikelihood', 1);
ll=0;
for i=1:model.numModels
    ll = ll + vargplvmLogLikelihood(model.comp{i});
end
f = (ll + KLdiv);



% Here we take advantage of the fact that in the bound, the likelihood and
% the KL part break into ll+KL. So, in the shared case, we only have
% KL + ll1 + ll2 + ...
function f = f2(model)

% We want to count for the KL part only once, because it is shared. Then, we
% will add all likelihood parts from all sub-models. The above can be
% achieved by calculating the bound for the first sub-model as usual, and
% then calculate only the likelihood part for the rest of the sub-models
% (because the KL part would be the same for all models since we called
% expandParam previously).
f = vargplvmLogLikelihood(model.comp{1}); % ll_1 + KL_1 (where KL_1 == KL_i for all i)
% The following call will not affect the model after returning from this
% function, since the model is not returned.
model = svargplvmPropagateField(model, 'onlyLikelihood', 1, true);
ll=0;
for i=2:model.numModels
    ll = ll + vargplvmLogLikelihood(model.comp{i}); % ll_i
end

f = f + ll;


