function model = vargpTimeDynamicsExpandParam(model, params)

% VARGPTIMEDYNAMICSEXPANDPARAM Place the parameters vector into the model for GP time dynamics.
%
%	Description:
%
%	MODEL = VARGPTIMEDYNAMICSEXPANDPARAM(MODEL, PARAMS) takes the given
%	vector of parameters and places them in the model structure, it then
%	updates any stored representations that are dependent on those
%	parameters, for example kernel matrices etc..
%	 Returns:
%	  MODEL - a returned model structure containing the updated
%	   parameters.
%	 Arguments:
%	  MODEL - the model structure for which parameters are to be
%	   updated.
%	  PARAMS - a vector of parameters for placing in the model
%	   structure.
%
%   COPYRIGHT : Andreas C. Damianou, 2010-2011
%   COPYRIGHT : Michalis K. Titsias, 2010-2011
%   COPYRIGHT : Neil D. Lawrence, 2010-2011
%
%	See also
%   vargpTimeDynamicsExtractParam, vargpTimeDynamicsCreate,
%   modelExpandParam
% VARGPLVM

% Parameters are received in the following order: (notation: % parameter{size})
% [dynamicsVardistParams{dynamics.vardist.nParams}        % mu_bar, lambda
%  dynamics.kernParams{dynamics.kern.nParams}]            % sigmaRBF, lt,  sigmaWhite


% variational parameters (means and covariances) reparametrized
startVal = 1;
endVal = model.vardist.nParams; 
model.vardist = modelExpandParam(model.vardist, params(startVal:endVal)); 
% dynamic model's kernel's parameters
startVal = endVal+1;
endVal = endVal + model.kern.nParams;
model.kern = kernExpandParam(model.kern, params(startVal:endVal));

    %%%% DEBUG_
    %fprintf(1,'In vargpTimeDynamicsExpandParam, params=%d %d %d\n', model.kern.comp{1}.inverseWidth, model.kern.comp{1}.variance,model.kern.comp{2}.variance );
    %%% _DEBUG

model.nParams = endVal;


