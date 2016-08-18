function [X, sigma2, W, model] = fgplvmEmbed(Y, dims, varargin)
% varargin is: {options, initVardistIters, iters, initSNR, display}
% e.g. Y = sin(1:0.05:10)' * rand(1,10); Y = Y + randn(size(Y)).*0.1;

iters = 200;
initSNR = 100;
display = 1;
if nargin > 3 && length(varargin)>5 && ~isempty(varargin{6})
    approx = varargin{6};
else
    approx = 'fitc';
end
options = fgplvmOptions(approx);
options.numActive = min(100, size(Y,1));
options.optimiser = 'scg';
    
if nargin > 3
    if ~isempty(varargin{1})
        options2 = varargin{1};
        options.numActive = options2.numActive;
        initSNR = options2.initSNR;
    end
    % initVardistIters is not used..
    if length(varargin)>1 && ~isempty(varargin{2}), initVardistIters = varargin{2}; end
    if length(varargin)>2 && ~isempty(varargin{3}), iters = varargin{3}; end
    if length(varargin)>3 && ~isempty(varargin{4}), initSNR = varargin{4}; end
    if length(varargin)>4 && ~isempty(varargin{5}), display = varargin{5}; end
end

latentDim = dims;
d = size(Y, 2);

model = fgplvmCreate(latentDim, d, Y, options);
%
if isfield(model, 'beta')
    model.beta = 1/((1/initSNR * var(model.m(:))));
    if model.beta < 1e-7
        warning('Beta was too small... Setting beta to 1e-7')
        model.beta = 1e-7;
    elseif model.beta > 1e+7
        warning('Beta was too big... Setting beta to 1e+7')
        model.beta = 1e+7;
    end
end


fprintf('#--- fgplvmEmbed from %d dims to %d dims, K=%d...\n',d,latentDim,options.numActive);

if initVardistIters > 0
    model.optimiseBeta = 0;
    fprintf('\n  # fgplvmEmbed: Initialising with fixed beta for %d iters...\n',initVardistIters);
    model = fgplvmOptimise(model, display, initVardistIters);
    model.optimiseBeta = 1;
end

if iters > 0
    fprintf('\n  # fgplvmEmbed: Optimising for %d iters...\n',iters);
    model = fgplvmOptimise(model, display, iters);
end
X = model.X;
sigma2 = [];
W = [];

if isfield(model, 'beta')
    finalSNR = var(model.m(:))*model.beta;
    fprintf('# fgplvmEmbed complete... final SNR: %f\n\n',finalSNR);
    if finalSNR < 10
        warning('# fgplvmEmbed: final SNR is too low!! (%f) \n\n', finalSNR)
    end
end

