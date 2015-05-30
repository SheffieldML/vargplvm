function [X, bias, scale] = outputsEmbed(Y, dims, scale2var1, scaleVal)
% OUTPUTSEMBED Initialise the latent space with the outputs themselves.
% VARGPLVM

if nargin < 4, scaleVal = [];   end
if nargin < 3, scale2var1 = 1; end
if nargin < 2 || isempty(dims)
    dims = size(Y,2);
elseif dims ~= size(Y,2)
    error('dims must equal output dimensionality.')
end
[X, bias, scale] = scaleData(Y, scale2var1, scaleVal);
