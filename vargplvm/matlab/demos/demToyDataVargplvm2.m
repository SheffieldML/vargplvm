randn('seed', 1e5);
rand('seed', 1e5);

N = 50;
D = 5;
if ~exist('mappingKern', 'var'), mappingKern = 'rbfardjit'; end
if ~exist('indPoints', 'var'), indPoints = 20; end
if ~exist('latentDim', 'var'), latentDim = 3; end
% Define a temporal model
if ~exist('dynamicsConstrainType', 'var'), dynamicsConstrainType = {'time'}; end
if ~exist('initVardistIters', 'var'), initVardistIters = 100; end
if ~exist('itNo', 'var'), itNo = [50 50]; end
if ~exist('dynamicKern', 'var'), dynamicKern = {'rbf','white','bias'}; end
if ~exist('makePlot', 'var'), makePlot = true; end

% Create Toy data
if ~exist('Y', 'var') && ~exist('t', 'var')
    t = linspace(0, 2*pi, N)';
    kernSamp = kernCreate(t, 'rbf');
    kernSamp.comp{1}.inverseWidth = 1;
    Ksamp = kernCompute(kernSamp, t);
    Y = real(gsamp(zeros(1,N), Ksamp, D))';
    Y = Y + randn(N,D)*0.05;
else
    N = size(Y,1);
    D = size(Y,2);
    assert(N == length(t));
end

%% Create and optimise the model

% Options for the model
vargplvm_init; % Returns a configuration structure 'globalOpt'
options = vargplvmOptions('dtcvar');
options.kern = mappingKern;
options.numActive = indPoints;
options.optimiser = 'scg2';
options.latentDim = latentDim;
options.initSNR = 100;

if ~isempty(dynamicsConstrainType)
    % Temporal model (VGPDS)
    if ~exist('optionsDyn', 'var')
        optionsDyn.type = 'vargpTime';
        optionsDyn.t=t;
        if ~isstruct(dynamicKern)
            optionsDyn.kern = kernCreate(t, dynamicKern);
        else
            optionsDyn.kern = dynamicKern;
        end
    end
    % Create and optimise the model
    [~, ~, ~, ~, model] = vargplvmEmbed2(Y, latentDim, options, initVardistIters, itNo, true, optionsDyn);
else
    % Non-dynamical model (Bayesian GP-LVM)
    % Create and optimise the model
    [~, ~, ~, ~, model] = vargplvmEmbed2(Y, latentDim, options, initVardistIters, itNo, true);
end

% Uncomment the following for to optimise for some more iterations
% model = vargplvmOptimiseModel(model, true, true, {initVardistIters, itNo}, true);

%% Plot the fit
if ~isempty(dynamicsConstrainType) && makePlot
    vgpdsPlotFit(model, [], 1:model.d, 6);
end
