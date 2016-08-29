% Augment the model with extra input-output pairs and optionally, inducing
% points.
function [model, ll]  = vargplvmAugmentModel(model, y_add, x_add, covs, X_u)

if nargin < 4 || isempty(covs)
    covs = 1e-11;
end
if nargin < 5, X_u = []; end


if isscalar(covs)
    covs = ones(size(x_add))*covs;
end

Nstar = size(y_add,1);

model.N = model.N + Nstar;
model.vardist.means = [model.vardist.means; x_add];
model.vardist.covars = [model.vardist.covars; covs];

if ~isempty(X_u)
    model.X_u = [model.X_u; X_u];
    model.k = model.k + size(X_u,1);
end


model.X = [model.X; x_add];
model.vardist.numData = model.N;
model.vardist.nParams = 2*prod(size(model.vardist.means));
model.vardist.transforms.index = model.N*model.q+1:model.vardist.nParams;
model.numParams = model.numParams + 2*model.q;
model.nParams = model.numParams;

% normalize y exactly as model.m is normalized
my = y_add - repmat(model.bias,size(y_add,1),1);
my = my./repmat(model.scale,size(y_add,1),1);

% change the data (by including the new point and taking only the present indices)
model.m = [model.m; my];
model.TrYY = sum(sum(model.m .* model.m));
if isfield(model, 'fixParamIndices')
    model.fixParamIndices = 1:2*model.N*model.q;
end
if isfield(model, 'fixInducing') && model.fixInducing
    model.inducingIndices = 1:model.k;
end


model = vargplvmUpdateStats(model, model.X_u);
if nargout > 1
    ll = vargplvmLogLikelihood(model);
end

