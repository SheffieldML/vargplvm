function X = svargplvmOptimiseDimVar(model,X,dim,display,iters,gradcheck,optimiser)

% SGPLVMOPTIMISEDIMVAR Optimise subspace of latent location
% FORMAT
% DESC Takes a fgplvm model and finds the latent location
% minimizing the variance on the output space
% ARG model : fgplvm model
% ARG X : latent initialisation
% ARG dim : dimensions to optimise
% ARG display : display optimisation iterations
% ARG iters : maximum number of iterations
% ARG gradcheck : check gradients
% ARG optimiser : optimiser (default = 'scg')
% RETURN X : optimised latent location
%
% SEEALSO : sgplvmOptimiseDimVarSequence
%
% COPYRIGHT : Neil D. Lawrence, Carl Henrik Ek, 2008
% Copied from the SGPLVM toolbox.
%
% SHEFFIELDML

if(nargin<7)
  optimiser = 'scg';
  if(nargin<6)
    gradcheck = false;
    if(nargin<5)
      iters = 100;
      if(nargin<4)
	display = 0;
	if(nargin<3)
	  error('Too Few Arguments');
	end
      end
    end
  end
end

% set optimisation settings
options = optOptions;
options(14) = iters;
options(9) = gradcheck;
options(1) = display;

optFunc = str2func(optimiser);

Ainv = model.P1' * model.P1; % size: NxN
model.Ainv = Ainv;
model.alpha = Ainv*model.Psi1'*model.m; % size: 1xD


for(i = 1:1:size(X,1))
  X(i,dim) = optFunc('varDimObjective',X(i,dim),options,'varDimGradient',model,X(i,:),dim);
end

return;
