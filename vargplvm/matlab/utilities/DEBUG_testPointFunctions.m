%save 'boundOLD.mat' boundOLD
%save 'modelOLD.mat' model
% At this point we assume that 'model' is a model acquired by running some
% demo (not for Point functions, for the whole dataset). 
%For the hardcoded params we have here, we ran the
% demStickTestPoint with q=10, M=50, N=50, with some other you need to
% change the hardcoded '20' which corresponds to the number of variational
% parameters mu and S (10*2)

clear;
DEBUGdemStickTestPoint % set: q=10, M=50, N=50

Nstar = 5;
params = vargplvmExtractParam(model);
model = vargplvmExpandParam(model, params);
Fold = vargplvmLogLikelihood(model);
Gold = vargplvmLogLikeGradients(model);

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

numIndPresent = round(0.5*model.d);
permi = randperm(model.d);
indexPresent =  permi(1:numIndPresent);
    
%ystar(:,indexPresent) = NaN; % With this in comments we should get zero difference in the result
%Fnew1 = vargplvmPointLogLikelihood(model2, vardistx, ystar);
Fnew = vargplvmPointLogLikelihood(model2, vardistx, ystar);
Fnew - Fold % This should be very small
%Gnew = vargplvmPointLogLikeGradient(model2, vardistx, ystar);
Gnew = vargplvmPointLogLikeGradient(model2, vardistx, ystar);
Goldmeans = reshape(Gold(1:model.N*model.q), model.N, model.q);
Goldcovs = reshape(Gold(model.N*model.q+1:model.N*model.q*2), model.N, model.q);
reorder = [Nstar+1:N-Nstar, 1:Nstar];
Goldmeans = Goldmeans(reorder,:);
Goldcovs = Goldcovs(reorder,:);
Gold1 = [Goldmeans(:)' Goldcovs(:)'];
max(abs(Gnew - Gold1))

vardistx = vardistCreate(model.vardist.means(1:Nstar,:), model.q, 'gaussian');
vardistx.means = model.dynamics.vardist.means(1:Nstar,:);
vardistx.covars = model.dynamics.vardist.covars(1:Nstar,:);

model2.vardistx = vardistx;
[x, varx, model2New] = vargplvmOptimisePoint(model2, vardistx, ystar, 1, 10);


