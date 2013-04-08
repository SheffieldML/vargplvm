%save 'boundOLD.mat' boundOLD
%save 'modelOLD.mat' model
Fold = vargplvmLogLikelihood(model);
model2 = model;
model2.y = model.y(2:end, :);
model2.m = model.m(2:end,:);
ystar = model.y(1,:);
model2.N = model.N - 1;
model2.dynamics.vardist.means = model.dynamics.vardist.means(2:end,:);
model2.dynamics.vardist.covars = model.dynamics.vardist.covars(2:end,:);
model2.dynamics.vardist.nParams = model.dynamics.vardist.nParams - 20;
model2.dynamics.t = model.dynamics.t(2:end);
vardistx.means = model.dynamics.vardist.means(1,:);
vardistx.covars = model.dynamics.vardist.covars(1,:);
vardistx.nParams = 20;
model2.vardist.means = model.vardist.means(2:end, : );
model2.vardist.covars = model.vardist.covars(2:end, : );
model2.vardist.nParams = 20;

model2.dynamics.t_star = model.dynamics.t(1);
model2.dynamics.N = model.dynamics.N-1;
[Fnew modelTemp] = vargplvmPointLogLikelihood(model2, vardistx, ystar);