function [X, sigma2, W, model, modelInitVardist] = vargplvmEmbed(Y, dims, varargin)
% VARGPLVMEMBED Embed given data Y into a lower dimensional structure with Bayesian GPLVM
% 
% COPYRIGHT: Andreas C. Damianou, 2012
% 
% VARGPLVM

% [X, sigma2, W, model, modelInitVardist] = vargplvmEmbed(Y, dims, varargin)
%
% varargin = {options, initIters, iters, display, optionsDyn}
% To do dynamical embedding, just call the function with all arguments in
% varargin (even if empty values are given, ie
% [..]=vargplvmEmbed(Y, dims, [],[],[],[],[]);

vargplvm_init;

initVardistIters = 15;
iters = 50;
display = 1;

options = vargplvmOptions('dtcvar');
options.kern = 'rbfardjit';
options.numActive = min(50,size(Y,1));
options.optimiser = 'scg2';
options.initSNR = 100;
dynUsed = 0;
optionsDyn = [];

if nargin > 2
    if ~isempty(varargin{1})
        options = varargin{1};
    end
    if length(varargin)>1 && ~isempty(varargin{2}), initVardistIters = varargin{2}; end
    if length(varargin)>2 && ~isempty(varargin{3}), iters = varargin{3}; end
    if length(varargin)>3 && ~isempty(varargin{4}), display = varargin{4}; end
    if length(varargin)>4 
        dynUsed=1;
        if ~isempty(varargin{5})
            optionsDyn = varargin{5}; 
        end
    end
end

if ~isfield(options, 'initSNR') && isfield(options, 'beta')
    globalOpt.betaInit = options.beta;
end

latentDim = dims;
d = size(Y, 2);

% demo using the variational inference method for the gplvm model
model = vargplvmCreate(latentDim, d, Y, options);
%
model = vargplvmParamInit(model, model.m, model.X); 
model = vargplvmModelInit(model, globalOpt);

fprintf('#--- vargplvmEmbed from %d dims to %d dims (initSNR=%f)...\n',d,latentDim,vargplvmShowSNR(model, false));


if dynUsed   
    %-------- Add dynamics to the model -----
    fprintf('  # Adding dynamics to the model...\n')
    if isempty(optionsDyn)
        optionsDyn.type = 'vargpTime';
        optionsDyn.t=[];
        optionsDyn.inverseWidth=30;
        optionsDyn.initX = model.X;%options.initX;
    end
    
    % Fill in with default values whatever is not already set
    optionsDyn = vargplvmOptionsDyn(optionsDyn, model.X);
    optionsDyn.initX = options.initX;
    model = vargplvmAddDynamics(model, 'vargpTime', optionsDyn, optionsDyn.t, 0, 0,optionsDyn.seq);
    model = vargplvmInitDynamics(model,optionsDyn);
end



modelInit = model;%%%% Delete

% Optimise the model.
fprintf('  # vargplvmEmbed: Optimising var. distr. for %d iters...\n',initVardistIters);
if initVardistIters > 0
    model.initVardist = 1; model.learnSigmaf = false;
    model = vargplvmOptimise(model, display, initVardistIters);
end
modelInitVardist = model;
model.initVardist = false; model.learnSigmaf = true;
if iters > 0
    fprintf('\n  # vargplvmEmbed: Optimising for %d iters...\n',iters);
    model = vargplvmOptimise(model, display, iters);
end
X = model.vardist.means;
sigma2 = model.vardist.covars;
W = vargplvmScales('get', model);

SNRfinal = vargplvmShowSNR(model, false);
fprintf('#--- Finished embedding. SNR: %f \n',SNRfinal);
fprintf('#--- Scales: %s\n\n', num2str(W));
if SNRfinal < 10
    warning(['During vargplvmEmbed SNR was too low (' num2str(SNRfinal) ')!'])
end
