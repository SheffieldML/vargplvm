function vardist = vardistExpandParam(vardist, params)

% VARDISTEXPANDPARAM Expand a parameter vector into a vardist structure.
% FORMAT
% DESC takes an VARDIST structure and a vector of parameters, and
% fills the structure with the given parameters. Also performs any
% necessary precomputation for likelihood and gradient
% computations, so can be computationally intensive to call.
% ARG model : the VARDIST structure to put the parameters in.
% ARG params : parameter vector containing the parameters to put in
% the VARDIST structure.
% 
% COPYRIGHT : Michalis K. Titsias, 2009
%
% COPYRIGHT : Neil D. Lawrence, 2009
%
% 
% SEEALSO : vardistCreate, vardistExtractParam, modelExpandParam

% VARGPLVM

if ~isempty(vardist.transforms)
  for i = 1:length(vardist.transforms)
    index = vardist.transforms(i).index;
    fhandle = str2func([vardist.transforms(i).type 'Transform']);
    params(index) = fhandle(params(index), 'atox');
  end
end

means = params(1:(vardist.numData*vardist.latentDimension));
st = vardist.numData*vardist.latentDimension + 1;
covs = params(st:end);

vardist.means = reshape(means, vardist.numData, vardist.latentDimension);
vardist.covars = reshape(covs, vardist.numData, vardist.latentDimension);