%DEMTOYDATAVARGPLVM1 Run Variational GPLVM with/without dynamics on toy data.
% DESC The toy data created are 3 trigonometric functions corrupted with
% noise, mapped to a 10-dimensional space and concatenated. The real
% dimensionality of this dataset is, thus, 3. This demo then runs
% variational GPLVM so that it manages to a) discover that the true
% dimensionality is 3 (this is obvious by looking at the final lengthscales)
% b) find latent functions that look like the original signal and
% c)can extrapolate in time and make predictions.
% It can be seen how the dynamics constrain the latent functions to be
% smooth (especially when the original signal is corrupted with a large
% amount of noise).
%
% Hints:
% You can experiment with the noiseLevel parameter in the
% vargplvmCrateToyData.m function, with different noise kernels for the dynamics,
% (i.e. 'white' or 'whitefixed', the later having a fixed variance), with/without
% the addition of dynamics, etc. Keep in mind that this is a trivial
% dataset, where the problem is actually linear and GPLVM is initialised
% with pca, so that -especially when no dynamics are used- the initial
% latent functions do not differ very much from the final ones.

% COPYRIGHT Andreas C. Damianou, 2011
% VARGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

% Constants
dataSetName = 'toy';
experimentNo = 1;
latentDim = 6; % Anything > 3 and < 10
dynamicsUsed = true;

vargplvmCreateToyData

% Set up model
options = vargplvmOptions('dtcvar');
options.kern = {'linard2', 'bias', 'white'};
options.numActive = 60; % Number of inducing points
options.optimiser = 'scg';
timeStampsTraining = t;

d = size(Y, 2);
fprintf(1,'# Creating the model...\n');
model = vargplvmCreate(latentDim, d, Y, options);
model = vargplvmParamInit(model, model.m, model.X);
model.beta=1/(0.01*var(model.m(:)));
model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));

if dynamicsUsed
    %-------- Add dynamics to the model -----
    optionsDyn.type = 'vargpTime';
    optionsDyn.t=t;

    % Dynamic kernel:
    % Putting a "whitefixed" instead of "white", will make the model give
    % samples x(t) that are more realisticly related to the observed data.
    kern = kernCreate(t, {'rbfperiodic', 'white'});
    % The following is related to the expected number of
    % zero-crossings.(larger inv.width numerator, rougher function)
    if ~strcmp(kern.comp{1}.type,'ou')
        kern.comp{1}.inverseWidth = 5./(((max(t)-min(t))).^2);
        kern.comp{1}.variance = 1;
    end
    optionsDyn.kern = kern;

    % Fill in with default values whatever is not already set
    optionsDyn = vargplvmOptionsDyn(optionsDyn);
    model = vargplvmAddDynamics(model, 'vargpTime', optionsDyn, optionsDyn.t, 0, 0,optionsDyn.seq);
    fprintf(1,'# Further calibration of the initial parameters...\n');
    model = vargplvmInitDynamics(model,optionsDyn);
end
modelInit = model;

%---------- OPTIMISATION ------
% do not learn beta for few iterations for intitilization
model.learnBeta = 0;
display = 1;
fprintf(1,'# Intitiliazing the model (fixed beta) %d iterations...\n',300);
model = vargplvmOptimise(model, display, 300);
disp('# Saving model after optimising beta...')
modelWriteResult(model, dataSetName, experimentNo);

% Optimise the model.
model.learnBeta = 1;
iters = 1000; % Default: 1000
fprintf(1,'# Optimising the model for %d iterations...\n',iters);
model = vargplvmOptimise(model, display, iters);
% Save the results.
fprintf(1,'# Saving model after doing %d iterations\n',iters)
modelWriteResult(model, dataSetName, experimentNo);


%--------- SIMPLE EVALUATION -----
% See the final lengthscales (see how some dimensions are switched-off).
bar(model.kern.comp{1}.inputScales); title('final lengthscales'); xlabel('q')

[a,ind]=sort(model.kern.comp{1}.inputScales,'descend'); % sort lengthscales
fprintf('Plotting the 3 latent functions corresponding to the largest 3 dimensions (scales):\n');
% Note: this is a plot of the latent functions x_q, they don't necessarily have
% to match the observed data Y, but their shape can provide insights, e.g. we
% know to expect something in the form of sines and cosines.
figure
for i=1:3
    subplot(1,3,i)
    plot(model.X(:,ind(i)));hold on, plot(Z{i},'r'), hold off
end
xlabel('N'); legend('latent function','original signal');

if dynamicsUsed
    figure
    dt = t(end)-t(end-1);
    timeStampsTest = (t(end)+dt:dt:t(end)+100*dt)'; % 100 time points in the future
    % Predict only given a test time vector (no partial information)
    [Testmeans2 Testcovars2] = vargplvmPredictPoint(model.dynamics, timeStampsTest);
    Varmu2 = vargplvmPosteriorMeanVar(model, Testmeans2, Testcovars2); % Predicted values
    model.dynamics.t_star = timeStampsTest;
    % Sample
    [ySamp, xSamp] = vargpTimeDynamicsSample(model, 1);
    d=10;
    plot(t,Y(:,d)), hold on, plot(timeStampsTest,Varmu2(:,d),'r'), plot(timeStampsTest,ySamp(:,d),'g')
    title(['Plot for dimension d=' num2str(d)]); legend('Y(:,d)','Prediction given t','Sampling'); xlabel('t');
end