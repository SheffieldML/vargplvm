%save 'boundOLD.mat' boundOLD
%save 'modelOLD.mat' model
% At this point we assume that 'model' is a model acquired by running some
% demo (not for Point functions, for the whole dataset). 
%For the hardcoded params we have here, we ran the
% demStickTestPoint with q=10, M=50, N=50, with some other you need to
% change the hardcoded '20' which corresponds to the number of variational
% parameters mu and S (10*2)

clear;

% Run with all blocks
DEBUGdemCmu35gplvmVargplvm1 % set: q=5, M=10, N=296 seq=[89 190 296]
seq=model.dynamics.seq;

% number of block out
Nblockout = 1; 

Nstar = seq(Nblockout);
params = vargplvmExtractParam(model);
model = vargplvmExpandParam(model, params);
Fold = vargplvmLogLikelihood(model);
Gold = vargplvmLogLikeGradients(model);

% Remove one block
model2 = model;
model2.y = model.y(Nstar+1:end, :);
model2.m = model.m(Nstar+1:end,:);
ystar = model.y(1:Nstar,:);
model2.N = model.N - Nstar;
model2.dynamics.vardist.means = model.dynamics.vardist.means(Nstar+1:end,:);
model2.dynamics.vardist.covars = model.dynamics.vardist.covars(Nstar+1:end,:);
model2.dynamics.vardist.nParams = model.dynamics.vardist.nParams - (2*model.q*Nstar);
model2.dynamics.t = model.dynamics.t(Nstar+1:end);
vardistx.means = model.dynamics.vardist.means(1:Nstar,:);
vardistx.covars = model.dynamics.vardist.covars(1:Nstar,:);
vardistx.nParams = 2*model.q*Nstar;
model2.vardist.means = model.vardist.means(Nstar+1:end, : );
model2.vardist.covars = model.vardist.covars(Nstar+1:end, : );
model2.vardist.nParams = size(model.vardist.means(Nstar+1:end,:),1);

model2.dynamics.t_star = model.dynamics.t(1:Nstar);
model2.dynamics.N = model.dynamics.N-Nstar;
model2.dynamics.seq=model2.dynamics.seq(Nblockout+1:end) - model2.dynamics.seq(Nblockout);



model2.dynamics.vardist.numData = model2.dynamics.N;
model2.vardist.numData = model2.N;
model2.vardist.transforms.type = model.vardist.transforms.type;
model2.vardist.transforms.index = (size(model2.vardist.means,2)*model2.N+1):size(model2.vardist.means,2)*model2.N*2;
model2.dynamics.vardist.transforms.index = (size(model2.vardist.means,2)*model2.N+1):size(model2.vardist.means,2)*model2.N*2;
model2.dynamics.nParams = model2.dynamics.kern.nParams + model2.dynamics.vardist.nParams;
%model2.dynamics.seq = model2.dynamics.N;
params = vargplvmExtractParam(model2);
model2.nParams = size(params,2);
model2 = vargplvmExpandParam(model2, params);

YtsOriginal = ystar;

numdataObs = round(0.5*Nstar);
permi = randperm(Nstar);
indexMissingData = permi(1:numdataObs);
numIndPresent = round(0.5*model.d);
permi = randperm(model.d);
indexPresent =  permi(1:numIndPresent);
indexMissing = setdiff(1:model.d, indexPresent);
    
ytmp = ystar(indexMissingData,:);
ytmp(:, indexMissing) = NaN;
ystar(indexMissingData,:) = ytmp;
%ystar(:, indexPresent) = NaN;
%Fnew = vargplvmPointLogLikelihood(model2, vardistx, ystar);
[Fnew, mm] = vargplvmSeqDynLogLikelihood(model2, vardistx, ystar);
fprintf(1,'# Fnew-Fold=%d \n',Fnew - Fold); % This should be very small
%Gnew = vargplvmPointLogLikeGradient(model2, vardistx, ystar);
Gnew = vargplvmSeqDynLogLikeGradient(model2, vardistx, ystar);
Goldmeans = reshape(Gold(1:model.N*model.q), model.N, model.q);
Goldcovs = reshape(Gold(model.N*model.q+1:model.N*model.q*2), model.N, model.q);
N=model.N; %%%???????
%reorder = [Nstar+1:N-Nstar, 1:Nstar];
Goldmeans = Goldmeans(1:Nstar,:);
Goldcovs = Goldcovs(1:Nstar,:);
Gold1 = [Goldmeans(:)' Goldcovs(:)'];
fprintf(1,'# Gnew-Gold=%d\n',max(abs(Gnew - Gold1)));


vardistx = vardistCreate(model.vardist.means(1:Nstar,:), model.q, 'gaussian');
vardistx.means = model.dynamics.vardist.means(1:Nstar,:);
vardistx.covars = model.dynamics.vardist.covars(1:Nstar,:);

model2.vardistx = vardistx;
%vargplvmOptimisePoint(model2, vardistx, ystar, 1, 1);
vargplvmOptimiseSeqDyn(model2, vardistx, ystar, 1, 10);

doPred = 1;
iters = 200;
if doPred   
  % Do also reconstruction in test data 
  [x, varx] = vargplvmOptimiseSeqDyn(model2, vardistx, ystar, 1, iters);
  % keep the optimized variational parameters
  barmu = x;
  lambda = varx;
 
  % Get the variational means and variances for the new test sequences and 
  % update the model to be prepared for prediction  
  [x, varx, modelUpdated] = vargplvmDynamicsUpdateModelTestVar(model2, barmu, lambda, ystar);
  
  % Latent variables corresponding to the data  with missing dimensions
  Testmeans = x(indexMissingData, :);
  Testcovars = varx(indexMissingData, :);
 
  YtsOrMs = YtsOriginal(indexMissingData,:);
  
  % Reconstruct all dimensions (missing and observed) for the test data  
  [mu, sigma] = vargplvmPosteriorMeanVar(modelUpdated, Testmeans, Testcovars);
  Varmu = mu; 
  Varsigma = sigma; 
  Varerrsum = sum((Varmu(:,indexMissing) - YtsOrMs(:,indexMissing)).^2);
  Varerror = mean(Varerrsum);
  
  % Error without updating the model structure (useful for de-bug... normally it must be less than 
  % Varerror)
  [muNoUpdated, sigmaNoUpdated] = vargplvmPosteriorMeanVar(model2, Testmeans, Testcovars);
  VarmuNoUpdated = muNoUpdated; 
  VarsigmaNoUpdated = sigmaNoUpdated; 
  VarerrsumNoUpdated = sum((VarmuNoUpdated(:,indexMissing) - YtsOrMs(:,indexMissing)).^2);
  VarerrorNoUpdated = mean(VarerrsumNoUpdated);
  
  % Random prediction 
  Randmu = sqrt(var(model2.m(:)))*randn(size(YtsOrMs(:,indexMissing)));
  % Tranform these predictions to the  initial space of y
  Randmu = Randmu.*repmat(model2.scale(indexMissing), size(Randmu,1),1);  
  Randmu = Randmu + repmat(model2.bias(indexMissing), size(Randmu,1),1); 
  Randerrsum = sum((Randmu - YtsOrMs(:,indexMissing)).^2);
  Randerror = mean(Randerrsum);
  
  % If you want to get the errors in TAYLOR's space 
  % then you need to know taylorbias and taylorscale. 
  % Then, you transfrom YtsOrMs, Varmu, Randmu. 
  %
  % For example:
  %
  % YtsOrMs = YtsOrMs - repmat(taylorbias, size(YtsOrMs, 1), 1);
  % YtsOrMs = YtsOrMs.*repmat(taylorscale, size(YtsOrMs, 1), 1);
  %
  % and similarly for Varmu etc. Then you compute the errors   
end
