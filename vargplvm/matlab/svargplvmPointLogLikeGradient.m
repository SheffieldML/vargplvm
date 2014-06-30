function [g, model] = svargplvmPointLogLikeGradient(model, vardistx, y)
% SVARGPLVMPOINTLOGLIKEGRADIENT Log-likelihood gradient for of some points of the MRD.
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
% SEEALSO : svargplvmPointLogLikelihood, svargplvmOptimisePoint
%
% COPYRIGHT : Andreas C. Damianou, 2013
% VARGPLVM


if isfield(model, 'dynamics') && ~isempty(model.dynamics)
    dynUsed = 1;
else
    dynUsed = 0;
end


% Indices of missing dimension
for i=model.testModalities
    indexMissing{i} = model.comp{i}.testPrecomp.indexMissing;
    indexPresent{i} = setdiff(1:model.comp{i}.d, indexMissing{i});
end

if dynUsed
    error('Not implemented for dynamcs')
    % use a special function
    %[g, model] = dynPointLogLikeGradient(model, vardistx, y, indexPresent, indexMissing);
    % g = g(:)'; %%% RE-OPT-CODE-REM
end

NstasrQ = vardistx.numData*vardistx.latentDimension;
gVarcovs0 = zeros(1, NstasrQ);
gVarcovs1 = zeros(1, NstasrQ);
gVarcovs2 = zeros(1, NstasrQ);
gVarmeans0 = zeros(1, NstasrQ);
gVarmeans1 = zeros(1, NstasrQ);
gVarmeans2 = zeros(1, NstasrQ);

% normalize y exactly as model.m is normalized
%my = y - repmat(model.bias(indexPresent),size(y,1),1);
%my = my./repmat(model.scale(indexPresent),size(y,1),1);

for i=model.testModalities
    [gPsi0, gPsi1, gPsi2] = vargpCovGrads(model.comp{i}, vardistx, y{i}, indexPresent{i});
    
    [gKern1, gVarmeans1tmp, gVarcovs1tmp] = kernVardistPsi1Gradient(model.comp{i}.kern, vardistx, model.comp{i}.X_u, gPsi1');
    [gKern2, gVarmeans2tmp, gVarcovs2tmp] = kernVardistPsi2Gradient(model.comp{i}.kern, vardistx, model.comp{i}.X_u, gPsi2);
    [gKern0, gVarmeans0tmp, gVarcovs0tmp] = kernVardistPsi0Gradient(model.comp{i}.kern, vardistx, gPsi0);
    
    gVarcovs0 = gVarcovs0 + (gVarcovs0tmp(:).*vardistx.covars(:))';
    gVarcovs1 = gVarcovs1 + (gVarcovs1tmp(:).*vardistx.covars(:))';
    gVarcovs2 = gVarcovs2 + (gVarcovs2tmp(:).*vardistx.covars(:))';
    
    gVarmeans0 = gVarmeans0 + gVarmeans0tmp;
    gVarmeans1 = gVarmeans1 + gVarmeans1tmp;
    gVarmeans2 = gVarmeans2 + gVarmeans2tmp;
end

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
%model.m = model.m(:,indexPresent);
%model.m = [my; model.m];

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
if isfield(model, 'DgtN') && model.DgtN
    model.P = model.P1 * (model.Psi1' * model.testPrecomp.mReducedGrad);
else
    model.P = model.P1 * (model.Psi1' * [my; model.m(:, indexPresent)]);
end
model.TrPP = sum(sum(model.P .* model.P));
model.B = model.P1' * model.P;
P1TP1 = (model.P1' * model.P1);
Tb = (1/model.beta) * d * P1TP1;
Tb = Tb + (model.B * model.B');
model.T1 = d * model.invK_uu - Tb;

gPsi2 = (model.beta/2) * model.T1;

gPsi0 = -0.5 * model.beta * d;

gPsi1 = model.beta*(P1TP1*model.Psi1'* model.testPrecomp.mY);

