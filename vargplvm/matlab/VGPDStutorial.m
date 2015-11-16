% VGPDSTUTORIAL Demonstration of VGPDS.
%
% COPYRIGHT: Andreas C. Damianou, 2014,2015
%
% VARGPLVM


% Run a subset (or all) of the scenarios (1 to 8).
% Give no arguments to the function to run all scenarios.
function model = VGPDStutorial(scenarios)
if nargin < 1
    scenarios = [1 2 3 4 5 6 7 8];
end

if ~isempty(intersect(1, scenarios))
    fprintf('# Scenario 1: default settings (Toy data, RBF cov. function for dynamics and mapping).\nPress any key to start...\n'); pause
    model = scenario1;
    fprintf('\n\n End of scenario 1!\n\n\n');
end

if ~isempty(intersect(2, scenarios))
    fprintf('# Scenario 2: As the previous scenario, but now do NOT use dynamics. The 2D projection will hence be less smooth (as X is no longer constrained to be a timeseries)\nPress any key to start...\n');pause
    model = scenario2;
    fprintf('\n\n End of scenario 2!\n\n\n');
end

if ~isempty(intersect(3, scenarios))
    fprintf('# Scenario 3:  Like the first scenario, but use a linear mapping from the latent space to the data space (the dynamics still have a non-linear mapping).\nPress any key to start...\n');pause
    model = scenario3;
    fprintf('\n\n End of scenario 3!\n\n\n');
end

if ~isempty(intersect(4, scenarios))
    fprintf('# Scenario 4: Both kernels (dynamics and mapping) are linear. There is no source of non-linearity and the result will be a straight line.\nPress any key to start...');pause
    model = scenario4;
    fprintf('\n\n End of scenario 4!\n\n\n');
end

if ~isempty(intersect(5, scenarios))
    fprintf('# Scenario 5:  Try a rougher kernel for the dynamics. Observe the stronger diagonal structure of the dyn. kernel compared to scenario 1.\nPress any key to start...\n');pause
    model = scenario5;
    fprintf('\n\n End of scenario 5!\n\n\n');
end

if ~isempty(intersect(6, scenarios))
    fprintf('# Scenario 6: Model a periodic signal.\nPress any key to start...\n');pause
    model = scenario6;
    fprintf('\n\n End of scenario 6!\n\n\n');
end

if ~isempty(intersect(7, scenarios))
    fprintf('# Scenario 7: Model a pseudo-periodic signal (periodic + random drift).\nPress any key to start...\n');pause
    model = scenario7;
    fprintf('\n\n End of scenario 7!\n\n\n');
end

if ~isempty(intersect(8, scenarios))
    fprintf('# Scenario 8: Model independent sequences: points within the same sequence are correlated and points across different sequences are not. Code is optimized to handle this case. Observe the block-diagonal structure of the dynamical cov. matrix.\nPress any key to start...\n');pause
    model = scenario8;
    fprintf('\n\n End of scenario 8 and end of Demo!\n\n\n');
end
end


% Helper visualisation function
function modelPlots(model)
if isfield(model, 'dynamics') && ~isempty(model.dynamics)
    % Visualise dynamical covariance matrix
    figure; imagesc(model.dynamics.Kt); title('Dynamics covariance matrix K_x')
end
% Visualise the ARD weights (which can switch off dimensions in X if
% they're not needed, i.e. automatic dimensionality detection)
figure; vargplvmShowScales(model); title('ARD weights for the latent dimensions detecting the latent dimensionality.')
% You can also plot the latent space in 2D, by projecting the latent space
% onto the two dominant dimensions
mm = vargplvmReduceModel(model,2); % order wrt to the inputScales
figure; ax = axes; lvmScatterPlot(mm, [], ax); title('The two dominant latent dimensions (ie 2D projection)')
end



%% Run toy demo with default settings
function model = scenario1()
close all;
demToyDataVargplvm2 % Check this script
% Signal to noise ratio shows how much information the model explained as
% noise. Ideally, we want this ratio as large as possible (close to 1 means
% that everything was modelled as noise.. perhaps too few training
% iterations.)
vargplvmShowSNR(model)
modelPlots(model)
% You can extract all the parameters of the model and their names.
% Parameters are the parts that are optimised. They might seem many, but
% bear in mind that most of them (X_u and var*) are variational rather than
% model parameters.
[params, names] = modelExtractParam(model);
names
end

%% As the previous demo, but now do NOT use dynamics. The 2D projection
% will hence be less smooth (as X is no longer constrained to be a
% timeseries)
function model = scenario2()
close all;
dynamicsConstrainType = {}; % This says that p(X) is just a unit Gaussian N(0,I)
demToyDataVargplvm2
modelPlots(model)
end

%% Like the first demo, but use a linear mapping from the latent space to the data
% space (the dynamics still have a non-linear mapping)
function model = scenario3()
close all;
mappingKern = {'linard2', 'white', 'bias'};
demToyDataVargplvm2
modelPlots(model)
end

%% Both kernels (dynamics and mapping) are linear. There's no source of non-linearity and the
% result will be a straight line
function model = scenario4()
close all;
initVardistIters = 500; % More iterations for initialising the var. distr.
indPoints = 5; % In theory, a linear kernel only needs 2 inducing points
mappingKern = {'linard2', 'white', 'bias'};
dynamicKern = {'lin', 'white', 'bias'};
demToyDataVargplvm2
modelPlots(model)
end


%% Try a 'rougher' kernel for the dynamics
function model = scenario5()
close all;
dynamicKern = {'matern32','white','bias'};
demToyDataVargplvm2
modelPlots(model)
end


%% Try a periodic signal
function model = scenario6()
close all;
randn('seed', 1e5);
rand('seed', 1e5);
N = 80; D = 3; indPoints = 25;
t = linspace(0, 5*pi, N)';
Y = sin(t)*randn(1,D) + randn(N,D)*0.05;
kk ={'rbfperiodic','white','bias'};
dynamicKern = kernCreate(t, kk);
dynamicKern.comp{1}.period = 5*pi;
demToyDataVargplvm2
modelPlots(model)
end

%% Try a pseudo-periodic signal
function model = scenario7()
close all;
randn('seed', 1e4);
rand('seed', 1e4);
N = 80; D = 3; indPoints = 25;
t = linspace(0, 6*pi, N)';
Y1 = sin(t)*randn(1,D);
kernSamp = kernCreate(t, 'rbf');
kernSamp.comp{1}.inverseWidth = 1;
Ksamp = kernCompute(kernSamp, t);
Y2 = real(gsamp(zeros(1,N), Ksamp, D))';
Y = Y1 + 0.3*Y2;
Y = Y + randn(N,D)*0.02;
kk ={'rbfperiodic','white','bias','rbf'};
dynamicKern = kernCreate(t, kk);
dynamicKern.comp{1}.period = 6*pi;
demToyDataVargplvm2
modelPlots(model)
end


%% Model Independent Sequences
function model = scenario8()
close all;
randn('seed', 1e4);
rand('seed', 1e4);
N = 80; D = 3; indPoints = 25;
tt = linspace(1, 20, 20)';
t = [tt; tt; tt; tt];
seq = [20, 40, 60, 80];

kernSamp = kernCreate(tt, 'rbf');
kernSamp.inverseWidth = 0.05;
Y = real(gsamp(zeros(1,length(tt)), kernCompute(kernSamp, tt), D))';
Y = [Y; real(gsamp(zeros(1,length(tt)), kernCompute(kernSamp, tt), D))'];
kernSamp = kernCreate(tt, 'rbfperiodic');
kernSamp.inverseWidth = 0.05;
Y = [Y; real(gsamp(zeros(1,length(tt)), kernCompute(kernSamp, tt), D))'];
Y = [Y; real(gsamp(zeros(1,length(tt)), kernCompute(kernSamp, tt), D))'];
Y = Y + randn(N,D)*0.02;

optionsDyn.type = 'vargpTime';
optionsDyn.t=t;
optionsDyn.kern = kernCreate(t, {'rbf','white','bias'});
optionsDyn.seq = seq;
makePlot = false;
demToyDataVargplvm2
modelPlots(model)
end