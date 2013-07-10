function kernVardistTest(kernType, numData, numIn)

% KERNVARDISTTEST description.

  
% VARGPLVM

x = randn(numData, numIn);
x2 = randn(numData/2, numIn);
    
kern = kernCreate(x,kernType);

params = 0.2*randn(1,kern.nParams)./sqrt(randn(1,kern.nParams).^2);
kern = kernExpandParam(kern, params);

vardist = vardistCreate(x, numIn, 'gaussian');

params = randn(1,(numData*numIn*2));
%means = x';
params(1:numData*numIn) = x(:)';
vardist = vardistExpandParam(vardist, params);

Psi2 = kernVardistPsi2Compute(kern, vardist, x2);

K = kernVardistPsi1Compute(kern, vardist, x2);
covGrad = ones(size(K));
[gKern, gVarmeans, gVarcovars, gInd] = kernVardistPsi1Gradient(kern, vardist, x2, covGrad);
%
KPsi2 = kernVardistPsi2Compute(kern, vardist, x2);
covGradPsi2 = ones(size(KPsi2));
[gKernPsi2, gVarmeansPsi2, gVarcovarsPsi2, gIndPsi2] = kernVardistPsi2Gradient(kern, vardist, x2, covGradPsi2);
%
psi0 = kernVardistPsi0Compute(kern, vardist);
[gKernPsi0, gVarmeansPsi0, gVarcovarsPsi0] = kernVardistPsi0Gradient(kern, vardist, 1);


epsilon = 1e-6;
paramskern = kernExtractParam(kern);
paramsvar = vardistExtractParam(vardist);
%xx2 = x2';
params = [paramskern paramsvar x2(:)'];
origParams = params;
for i = 1:length(params);
  params = origParams;
  params(i) = origParams(i) + epsilon;
  kern = kernExpandParam(kern, params(1:kern.nParams));
  vardist = vardistExpandParam(vardist, params(kern.nParams+1:kern.nParams+vardist.nParams));
  xx2 = params(kern.nParams+vardist.nParams+1:end);
  x2 = reshape(xx2,[numData/2 numIn]);
  Lplus(i) = full(sum(sum(kernVardistPsi1Compute(kern, vardist, x2))));
  params(i) = origParams(i) - epsilon;
  kern = kernExpandParam(kern, params(1:kern.nParams));
  vardist = vardistExpandParam(vardist, params(kern.nParams+1:kern.nParams+vardist.nParams));
  xx2 = params(kern.nParams+vardist.nParams+1:end);
  x2 = reshape(xx2,[numData/2 numIn]);
  Lminus(i) = full(sum(sum(kernVardistPsi1Compute(kern, vardist, x2))));
end

params = origParams;
gLDiff = .5*(Lplus - Lminus)/epsilon;
g = [gKern gVarmeans gVarcovars gInd];

% check firstly the kernel hyperparameters 
kerndiff = abs(g(1:kern.nParams) - gLDiff(1:kern.nParams));
[g(1:kern.nParams); gLDiff(1:kern.nParams)]
pause
index = [kern.nParams+1:kern.nParams+(vardist.nParams/2)];
varmeansdiff =  abs(g(index) - gLDiff(index));
[g(index); gLDiff(index)]
pause
index = [kern.nParams+(vardist.nParams/2)+1:kern.nParams+vardist.nParams];
varcovarsdiff =  abs(g(index) - gLDiff(index));
[g(index); gLDiff(index)]
pause
index = [(kern.nParams+vardist.nParams+1):size(g,2)]; 
varinddiff =  abs(g(index) - gLDiff(index)); 
[g(index); gLDiff(index)]
pause
fprintf('Kernel hyperparameters\n');
if strcmp(kern.type,'rbfard2')
fprintf('var: %2.6g max and min inputscale: %2.6g, %2.6g.\n',kern.variance,max(kern.inputScales),min(kern.inputScales));
elseif strcmp(kern.type,'linard2')
fprintf('max and min inputscale: %2.6g, %2.6g.\n',max(kern.inputScales),min(kern.inputScales)); 
else
  % do nothing
end

fprintf('----- Psi1 term ------- \n');
fprintf('Kernel hyps max diff: %2.6g.\n', max(kerndiff));
fprintf('Variational means max diff: %2.6g.\n', max(varmeansdiff));
fprintf('Variational covars max diff: %2.6g.\n', max(varcovarsdiff));
fprintf('Inducing inputs max diff: %2.6g.\n', max(varinddiff));
fprintf('\n');


fprintf('----- Psi2 term ------- \n');
E = eig(Psi2); 
fprintf('max and min eigenvalue of Psi2: %2.6g, %2.6g.\n',max(E),min(E));
origParams = params;
for i = 1:length(params);
  params = origParams;
  params(i) = origParams(i) + epsilon;
  kern = kernExpandParam(kern, params(1:kern.nParams));
  vardist = vardistExpandParam(vardist, params(kern.nParams+1:kern.nParams+vardist.nParams));
  xx2 = params(kern.nParams+vardist.nParams+1:end);
  x2 = reshape(xx2,[numData/2 numIn]);
  Lplus(i) = full(sum(sum(kernVardistPsi2Compute(kern, vardist, x2))));
  params(i) = origParams(i) - epsilon;
  kern = kernExpandParam(kern, params(1:kern.nParams));
  vardist = vardistExpandParam(vardist, params(kern.nParams+1:kern.nParams+vardist.nParams));
  xx2 = params(kern.nParams+vardist.nParams+1:end);
  x2 = reshape(xx2,[numData/2 numIn]);
  Lminus(i) = full(sum(sum(kernVardistPsi2Compute(kern, vardist, x2))));
end
params = origParams;
gLDiff = .5*(Lplus - Lminus)/epsilon;
g = [gKernPsi2 gVarmeansPsi2 gVarcovarsPsi2 gIndPsi2];
% check firstly the kernel hyperparameters 
kerndiff = abs(g(1:kern.nParams) - gLDiff(1:kern.nParams));
[g(1:kern.nParams); gLDiff(1:kern.nParams)]
pause
index = [kern.nParams+1:kern.nParams+(vardist.nParams/2)];
varmeansdiff =  abs(g(index) - gLDiff(index));
[g(index); gLDiff(index)]
pause
index = [kern.nParams+(vardist.nParams/2)+1:kern.nParams+vardist.nParams];
varcovarsdiff =  abs(g(index) - gLDiff(index));
[g(index); gLDiff(index)]
pause
index = [(kern.nParams+vardist.nParams+1):size(g,2)]; 
varinddiff =  abs(g(index) - gLDiff(index)); 
[g(index); gLDiff(index)]
pause
fprintf('Kernel hyps max diff: %2.6g.\n', max(kerndiff));
fprintf('Variational means max diff: %2.6g.\n', max(varmeansdiff));
fprintf('Variational covars max diff: %2.6g.\n', max(varcovarsdiff));
fprintf('Inducing inputs max diff: %2.6g.\n', max(varinddiff));
fprintf('\n');

fprintf('----- Psi0 term ------- \n');
E = eig(Psi2); 
origParams = params;
for i = 1:length(params);
  params = origParams;
  params(i) = origParams(i) + epsilon;
  kern = kernExpandParam(kern, params(1:kern.nParams));
  vardist = vardistExpandParam(vardist, params(kern.nParams+1:kern.nParams+vardist.nParams));
  xx2 = params(kern.nParams+vardist.nParams+1:end);
  x2 = reshape(xx2,[numData/2 numIn]);
  Lplus(i) = full(sum(sum(kernVardistPsi0Compute(kern, vardist))));
  params(i) = origParams(i) - epsilon;
  kern = kernExpandParam(kern, params(1:kern.nParams));
  vardist = vardistExpandParam(vardist, params(kern.nParams+1:kern.nParams+vardist.nParams));
  xx2 = params(kern.nParams+vardist.nParams+1:end);
  x2 = reshape(xx2,[numData/2 numIn]);
  Lminus(i) = full(sum(sum(kernVardistPsi0Compute(kern, vardist))));
end
params = origParams;
gLDiff = .5*(Lplus - Lminus)/epsilon;
g = [gKernPsi0 gVarmeansPsi0 gVarcovarsPsi0];
% check firstly the kernel hyperparameters 
kerndiff = abs(g(1:kern.nParams) - gLDiff(1:kern.nParams));
[g(1:kern.nParams); gLDiff(1:kern.nParams)]
pause
index = [kern.nParams+1:kern.nParams+(vardist.nParams/2)];
varmeansdiff =  abs(g(index) - gLDiff(index));
[g(index); gLDiff(index)]
pause
index = [kern.nParams+(vardist.nParams/2)+1:kern.nParams+vardist.nParams];
varcovarsdiff =  abs(g(index) - gLDiff(index));
[g(index); gLDiff(index)]
pause
fprintf('Kernel hyps max diff: %2.6g.\n', max(kerndiff));
fprintf('Variational means max diff: %2.6g.\n', max(varmeansdiff));
fprintf('Variational covars max diff: %2.6g.\n', max(varcovarsdiff));

