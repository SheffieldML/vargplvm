function model = vargplvmExpandParam(model, params)

% VARGPLVMEXPANDPARAM Expand a parameter vector into a GP-LVM model.
% FORMAT
% DESC takes an VARGPLVM structure and a vector of parameters, and
% fills the structure with the given parameters. Also performs any
% necessary precomputation for likelihood and gradient
% computations, so can be computationally intensive to call.
% ARG model : the VARGPLVM structure to put the parameters in.
% ARG params : parameter vector containing the parameters to put in
% the VARGPLVM structure.
% 
%
% COPYRIGHT : Michalis K. Titsias, 2009-2011
% COPYRIGHT : Neil D. Lawrence, 2009-2011
% 
% Modifications: Andreas C. Damianou, 2010-2011 
%
% SEEALSO : vargplvmCreate, vargplvmExtractParam, modelExpandParam

% VARGPLVM


%%% Parameters must be passed as a vector in the following order (left to right) 
% - parameter{size} -
% vardistParams{model.vardist.nParams} % mu, S
%       OR
% [dynamicsVardistParams{dynamics.vardist.nParams} dynamics.kernParams{dynamics.kern.nParams}] % mu_bar, lambda
% inducingInputs{model.q*model.k}
% kernelParams{model.kern.nParams}
% beta{prod(size(model.beta))}


startVal = 1;
if isfield(model, 'dynamics') & ~isempty(model.dynamics)
    % variational parameters (reparametrized) AND dyn.kernel's parameters
    endVal = model.dynamics.nParams;
    model.dynamics = modelExpandParam(model.dynamics, params(startVal:endVal));
else
    % variational parameters (means and covariances), original ones
    endVal = model.vardist.nParams;
    model.vardist = modelExpandParam(model.vardist, params(startVal:endVal)); 
end

% inducing inputs: if fixInducing is active, they are tied to the
% variational means. If learnInducing exists as a field and is false, they
% are not included at all (they are just fixed parameters).
startVal = endVal+1;
if model.fixInducing  || (isfield(model, 'learnInducing') && ~model.learnInducing)
    if ~(isfield(model, 'learnInducing') && ~model.learnInducing) % If this is true, don't change the values at all
        if isfield(model, 'dynamics') & ~isempty(model.dynamics)
            Kt = kernCompute(model.dynamics.kern, model.dynamics.t);%%%%%%%%%%%%5
            model.X_u = Kt*model.dynamics.vardist.means; %dynamics
            model.X_u = model.X_u(model.inducingIndices,:);
        else
             model.X_u = model.vardist.means(model.inducingIndices, :); % static
        end
    end
    % X_u values are taken from X values.
   % model.X_u = model.X(model.inducingIndices, :);
else
    % Parameters include inducing variables.
    endVal = endVal + model.q*model.k;
    model.X_u = reshape(params(startVal:endVal),model.k,model.q);
end

% kernel hyperparameters 
startVal = endVal+1; 
endVal = endVal + model.kern.nParams;
model.kern = kernExpandParam(model.kern, params(startVal:endVal));

% likelihood beta parameters
if model.optimiseBeta
  startVal = endVal + 1;
  endVal = endVal + prod(size(model.beta));
  if ~isstruct(model.betaTransform)
    fhandle = str2func([model.betaTransform 'Transform']);
    model.beta = fhandle(params(startVal:endVal), 'atox');
  else
      if isfield(model.betaTransform,'transformsettings') && ~isempty(model.betaTransform.transformsettings)
          fhandle = str2func([model.betaTransform.type 'Transform']);
          model.beta = fhandle(params(startVal:endVal), 'atox', model.betaTransform.transformsettings);
      else
          error('vargplvmExtractParam: Invalid transform specified for beta.'); 
      end
  end
end

model.nParams = endVal;

% Update statistics
model = vargplvmUpdateStats(model, model.X_u);

% %%%TEMP: This is not needed, probably. If yes, it should be merged with
% %%%the above code for fixInducing.
% if model.fixInducing
%     model.X_u=model.X(model.inducingIndices, :);
% end
% %%%

