function model = vargpTimeDynamicsCreate(q, d, latentVals, options, varargin)

% VARGPTIMEDYNAMICSCREATE Create the time dynamics variational model.
%
%	Description:
%
%	MODEL = VARGPTIMEDYNAMICSCREATE(Q, Q, X, OPTIONS, T, DIFF, LEARN, SEQ)
%	creates a Gaussian process model for dealing with dynamics in the
%	latent space of a GP-LVM. The input to the dynamics model is time,
%	rather than the previous observation as is the case for GPDYNAMICS.
%	 Returns:
%	  MODEL - model structure containing the Gaussian process for the Dynamics.
%	 Arguments:
%	  Q - the latent space dimension.
%	  Q - the latent space dimension.
%	  X - the latent variables.
%	  OPTIONS - options structure as defined by gpOptions.m.
%	  T - the input time values.
%	  DIFF - Whether or not to use differences between points in the
%	   latent space as the targets for the GP or absolute location of
%	   points (default 1). --not implemented yet
%	  LEARN - Whether or not to learn the parameters of the dynamics
%	   model (default 0). -- not implemented yet
%	  SEQ - array containing the indices of the last frame of each
%	   sequence in the data. The array should be the same length as the
%	   number of sequences in the data (default []).
%
%   COPYRIGHT : Andreas C. Damianou, 2010-2011
%   COPYRIGHT : Michalis K. Titsias, 2010-2011
%   COPYRIGHT : Neil D. Lawrence, 2010-2011
%
%	See also
%	vargplvmCreate, vargplvmDynamicsCreate, vargpTimeDynamicsLogLikelihood


%	Based on GPTIMEDYNAMICSCREATE
% VARGPLVM

% Varargin may contain t, diff, learn, seq).
if nargin>4
  t = varargin{1};
else
  t = [1:size(latentVals, 1)]'; % The upper levels ensure that t will always be set
end
if nargin>5 % Not implemented yet
  diff = varargin{2};
else 
  diff = 0;
end
if nargin > 6 % Not implemented yet
  learn = varargin{3};
else
  learn = 0;
end
if nargin > 7 % Not implemented yet
  seq = varargin{4};
else
  %seq = size(latentVals, 1);
  seq=[];
end
  
if(iscell(options))
  varargin = options;
  options = varargin{1};
  t = varargin{2};
  
  % Not implemented yet
  diff = varargin{3}; 
  learn = varargin{4};
  seq = varargin{5};
end

if ~isfield(options, 'isMissingData')
    model.isMissingData = 0;
else
    model.isMissingData = options.isMissingData;
end    


if d ~= q
    error(['problems with q not equal to d cannot be solved with current implementation']);
end

model.t = t;
% This is (?) needed for visualisation functions
model.X = latentVals;


model.q = size(latentVals, 2);
% model.d = model.q;
model.N = size(latentVals, 1);

% if isfield(options, 'whiteVariance')
%     model.whiteVariance = options.whiteVariance;
% end



if size(t) ~= model.N
    error(['Input vector t does not have dimension ' num2str(model.N)]);
end



if any(isnan(t)) & ~options.isMissingData
  error('NaN values in t.')
end


if isstruct(options.kern) 
  model.kern = options.kern;
else
  model.kern = kernCreate(t, options.kern); 
end


% Initialize the model
if isfield(options, 'notransform') && options.notransform == true
    model.vardist = vardistCreate(latentVals, q, 'gaussian', 'identity');
else
    model.vardist = vardistCreate(latentVals, q, 'gaussian');
end
    

% Must find a formula which, given Kt, relates initial lambdas to initial Sqs
% so that Sqs are around 0.5
if ~isfield(options, 'lambda_init') 
    model.vardist.covars = 0.25 * ones(size(model.vardist.covars));
else
    model.vardist.covars = options.lambda_init * ones(size(model.vardist.covars));
end

model.vardist.means = latentVals; % vardistCreate does that

model.type = 'vargpTimeDynamics';
model.dynamicsType = 'regressive';

initParams = vargpTimeDynamicsExtractParam(model); 
% This forces kernel computation.
model = vargpTimeDynamicsExpandParam(model, initParams);


 model.seq = seq;
 if max(seq) > model.N
    error('The maximum seq. index cannot exceed N.')
 end
 if sum(seq>0) ~= length(seq)
     error('Sequence indices must be positive.')
 end
 
 % % Not implemented yet
% model.diff = diff;
% model.learn = learn;


