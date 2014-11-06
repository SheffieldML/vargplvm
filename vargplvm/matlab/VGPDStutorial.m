% VGPDSTUTORIAL Demonstration of VGPDS.
%
% COPYRIGHT: Andreas C. Damianou, 2014
%
% VARGPLVM

%-----------------    SIMPLE INITIAL DEMO -------------------------------%

%% Run toy demo with default settings
clear
demToyDataVargplvm2 % Check this script
% Visualise dynamical covariance matrix
figure
imagesc(model.dynamics.Kt); title('K_x')
% Visualise the ARD weights (which can switch off dimensions in X if
% they're not needed, i.e. automatic dimensionality detection)
figure
vargplvmShowScales(model);
% Signal to noise ratio shows how much information the model explained as
% noise. Ideally, we want this ratio as large as possible (close to 1 means
% that everything was modelled as noise.. perhaps too few training
% iterations.)
vargplvmShowSNR(model)
figure; ax = axes;
% You can also plot the latent space in 2D, by projecting the latent space
% onto the two dominant dimensions
mm = vargplvmReduceModel(model,2); % order wrt to the inputScales 
lvmScatterPlot(mm, [], ax); title('2D projection')
% You can extract all the parameters of the model and their names.
% Parameters are the parts that are optimised. They might seem many, but
% bear in mind that most of them (X_u and var*) are variational rather than
% model parameters.
[params, names] = modelExtractParam(model);

%% As the previous demo, but now do NOT use dynamics. The 2D projection
% will hence be less smooth (as X is no longer constrained to be a
% timeseries)
clear
dynamicsConstrainType = {}; % This says that p(X) is just a Gaussian N(0,I)
demToyDataVargplvm2
vargplvmShowScales(model);
figure; ax = axes;
mm = vargplvmReduceModel(model,2); % order wrt to the inputScales 
lvmScatterPlot(mm, [], ax); title('2D projection')

%% Like the first demo, but use a linear mapping from the latent space to the data
% space (the dynamics still have a non-linear mapping)
clear
mappingKern = {'linard2', 'white', 'bias'};
demToyDataVargplvm2

%% Both kernels are linear. There's no source of non-linearity and the
% result will be a straight line
clear
initVardistIters = 500; % More iterations for initialising the var. distr.
indPoints = 5; % In theory, a linear kernel only needs 2 inducing points
mappingKern = {'linard2', 'white', 'bias'};
dynamicKern = {'lin', 'white', 'bias'};
demToyDataVargplvm2

%% Try a 'rougher' kernel for the dynamics
clear
dynamicKern = {'matern32','white','bias'};
demToyDataVargplvm2

%% Try a periodic signal
clear
randn('seed', 1e5);
rand('seed', 1e5);
N = 80; D = 3; indPoints = 25;
t = linspace(0, 5*pi, N)';
Y = sin(t)*randn(1,D) + randn(N,D)*0.05;
kk ={'rbfperiodic','white','bias'};
dynamicKern = kernCreate(t, kk);
dynamicKern.comp{1}.period = 5*pi;
demToyDataVargplvm2


%% Try a pseudo-periodic signal
clear
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



%% Model Independent Sequences
clear
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

imagesc(model.dynamics.Kt); title('K_x')
