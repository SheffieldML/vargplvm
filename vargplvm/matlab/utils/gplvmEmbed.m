function [X, model] = gplvmEmbed(Y, dims, iters, options, numActive, display)
% GPLVMEMBED Embed given data Y into a lower dimensional structure with Bayesian GPLVM
% 
% COPYRIGHT: Andreas C. Damianou, 2013
% 
% VARGPLVM

% [X, model] = gplvmEmbed(Y, dims, iters, options ,numActive, display)


if nargin < 2 || isempty(dims), dims = 2; end
if nargin < 3 || isempty(iters), iters = 100; end
if nargin < 4 || isempty(options), options = fgplvmOptions('fitc'); end
if nargin < 5 || isempty(numActive), numActive = min(50, size(Y,1)); end
if nargin < 6 || isempty(display), display = 1; end

options.numActive = numActive;
d = size(Y, 2);
fprintf('#--- gplvmEmbed from %d dims to %d dims...\n',d,dims);
options.optimiser = 'scg2';
model = fgplvmCreate(dims, d, Y, options);
model = fgplvmOptimise(model, display, iters);
X = model.X;