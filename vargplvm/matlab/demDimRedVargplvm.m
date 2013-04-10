% DEMDIMREDVARGPLVM A simple demonstration of dimensionality reduction for
% the Bayesian GP-LVM
% COPYRIGHT: Andreas C. Damianou, 2012
% VARGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

% Constants
dataSetName = 'toy';
experimentNo = 1;
latentDim = 2; % Anything > 3 and < 10
dynamicsUsed = false;
embedLinearly = false;

%%%%%%%%%%%%%%
t = linspace(0,4*pi,100);
% The original signal will be a cosine, a sine and a squared cosine.
Z1 = cos(t)';
Z2 = sin(t)';
%Z3= (cos(t)').^2;

% Store the original signals. Normally, the corresponding latent functions
% found for these signals will not be ordered, i.e. Z1 does not necessarily
% correspond to x_1(t). But with fixed random seeds it seems that the following
% correspondence is done Z1 ->x_3(t), Z2->x_2(t), X3 ->x_1(t)
%if dynamicsUsed
%    Z{1} = Z3; Z{2} = Z2; Z{3} = Z1;
%else
%    Z{1} = Z1; Z{2} = Z2; Z{3} = Z3;
%end
Z{1} = Z1; Z{2} = Z2;

% Scale and center data
bias_Z1 = mean(Z1);
Z1 = Z1 - repmat(bias_Z1,size(Z1,1),1);
scale_Z1 = max(max(abs(Z1)));
Z1 = Z1 ./scale_Z1;

bias_Z2 = mean(Z2);
Z2 = Z2 - repmat(bias_Z2,size(Z2,1),1);
scale_Z2 = max(max(abs(Z2)));
Z2 = Z2 ./ scale_Z2;



noiseLevel = 0; % Default: 0.1 (or 0.5)


%------- LINEAR EMBEDING-----
% Map 1-D to 3-D and add some noise
if embedLinearly
    Z2p = Z2*rand(1,3);
    Z2p = Z2p + noiseLevel.*randn(size(Z2p));
    Z1p = Z1*rand(1,3);
    Z1p = Z1p + noiseLevel.*randn(size(Z1p));
else
    Z2p = exp((Z2*rand(1,3))).^2;
    Z2p = Z2p + noiseLevel.*randn(size(Z2p));
    Z1p = exp((Z1*rand(1,3))).^2;
    Z1p = Z1p + noiseLevel.*randn(size(Z1p));
%    x = linspace(-1, 1, 100)'; % input indices in the x-axis
%     trueKern = kernCreate(x, 'matern32'); % Generate samples from this kernel
%     K = kernCompute(trueKern, x) + eye(size(x, 1))*noiseVar;
%     % Sample some true function values.
%     yTrue = gsamp(zeros(size(x))', K, 3)';
end


% Form dataset by concatenating the signals
Y = [Z1p Z2p];
t = t';
%%%%%%%%%%%%%%

% Set up model
options = vargplvmOptions('dtcvar');
if embedLinearly
    options.kern = {'linard2', 'bias', 'white'};
else
    options.kern = {'rbfard2', 'bias', 'white'};
end
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
%modelWriteResult(model, dataSetName, experimentNo);

% Optimise the model.
model.learnBeta = 1;
iters = 200; % Default: 1000
fprintf(1,'# Optimising the model for %d iterations...\n',iters);
model = vargplvmOptimise(model, display, iters);
% Save the results.
fprintf(1,'# Saving model after doing %d iterations\n',iters)
%modelWriteResult(model, dataSetName, experimentNo);


%--------- SIMPLE EVALUATION -----
% See the final lengthscales (see how some dimensions are switched-off).
bar(model.kern.comp{1}.inputScales); title('final lengthscales'); xlabel('q')

[a,ind]=sort(model.kern.comp{1}.inputScales,'descend'); % sort lengthscales
fprintf('Plotting the 2 latent functions corresponding to the largest 3 dimensions (scales):\n');
% Note: this is a plot of the latent functions x_q, they don't necessarily have
% to match the observed data Y, but their shape can provide insights, e.g. we
% know to expect something in the form of sines and cosines.
figure
for i=1:2
    subplot(1,2,i)
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
    d=3;
    plot(t,Y(:,d)), hold on, plot(timeStampsTest,Varmu2(:,d),'r'), plot(timeStampsTest,ySamp(:,d),'g')
    title(['Plot for dimension d=' num2str(d)]); legend('Y(:,d)','Prediction given t','Sampling'); xlabel('t');
end



%%%%%%%%%%%%%%%%%%%%%
%pcaX = pcaEmbed(Y,2);
[U,V] = pca(Y,2);
pcaX = Y*V;

subplot(4,2,1); plot3(model.X(:,1), model.X(:,2),t); xlabel('x_1'); ylabel('x_2'); zlabel('t');
subplot(4,2,2); plot3(Z1,Z2,t);     xlabel('Z1'); ylabel('Z2'); zlabel('t');
subplot(4,2,3); plot(model.X(:,1)); ylabel('x_1'); xlabel('t');
subplot(4,2,4); plot(model.X(:,2)); ylabel('x_2'); xlabel('t');
subplot(4,2,5); plot(pcaX(:,1)); ylabel('pcaX_1'); xlabel('t');
subplot(4,2,6); plot(pcaX(:,2)); ylabel('pcaX_2'); xlabel('t');
subplot(4,2,7); surf(Z1p);          title('Z1p');
subplot(4,2,8); surf(Z2p);          title('Z2p');

