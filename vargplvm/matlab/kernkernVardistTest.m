function kernkernVardistTest(kernType1, kernType2, numData, numIn)

% KERNKERNVARDISTTEST Description

% VARGPLVM

x = randn(numData, numIn);
x2 = randn(numData/2, numIn);
    
kern1 = kernCreate(x,kernType1);
kern2 = kernCreate(x,kernType2);

params1 = 0.2*randn(1,kern1.nParams)./sqrt(randn(1,kern1.nParams).^2);
kern1 = kernExpandParam(kern1, params1);

params2 = 0.2*randn(1,kern2.nParams)./sqrt(randn(1,kern2.nParams).^2);
kern2 = kernExpandParam(kern2, params2);


vardist = vardistCreate(x, numIn, 'gaussian');

params = randn(1,(numData*numIn*2));
params(1:numData*numIn) = x(:)';
vardist = vardistExpandParam(vardist, params);

Psi2 = kernkernVardistPsi2Compute(kern1, kern2, vardist, x2);
covGrad = ones(size(Psi2));
[gKern1, gKern2, gVarmeans, gVarcovars, gInd] = kernkernVardistPsi2Gradient(kern1, kern2, vardist, x2, covGrad);
%

epsilon = 1e-6;
paramskern1 = kernExtractParam(kern1);
paramskern2 = kernExtractParam(kern2);
paramsvar = vardistExtractParam(vardist);
%xx2 = x2';
params = [paramskern1 paramskern2 paramsvar x2(:)'];
origParams = params;

fprintf('Kernel hyperparameters\n');
if strcmp(kern1.type,'rbfard2')
fprintf('var: %2.6g max and min inputscale: %2.6g, %2.6g.\n',kern1.variance,max(kern1.inputScales),min(kern1.inputScales));
elseif strcmp(kern1.type,'linard2')
fprintf('max and min inputscale: %2.6g, %2.6g.\n',max(kern1.inputScales),min(kern1.inputScales)); 
else
  % do nothing
end

fprintf('----- Psi2 term ------- \n');
E = eig(Psi2); 
fprintf('max and min eigenvalue of Psi2: %2.6g, %2.6g.\n',max(E),min(E));
origParams = params;
for i = 1:length(params);
  params = origParams;
  params(i) = origParams(i) + epsilon;
  kern1 = kernExpandParam(kern1, params(1:kern1.nParams));
  kern2 = kernExpandParam(kern2, params(kern1.nParams+1:kern1.nParams+kern2.nParams));
  
  vardist = vardistExpandParam(vardist, params(kern1.nParams+kern2.nParams+1:kern1.nParams+kern2.nParams+vardist.nParams));
  xx2 = params(kern1.nParams+kern2.nParams+vardist.nParams+1:end);
  x2 = reshape(xx2,[numData/2 numIn]);
  Lplus(i) = full(sum(sum(kernkernVardistPsi2Compute(kern1, kern2, vardist, x2))));
  params(i) = origParams(i) - epsilon;
  
  kern1 = kernExpandParam(kern1, params(1:kern1.nParams));
  kern2 = kernExpandParam(kern2, params(kern1.nParams+1:kern1.nParams+kern2.nParams));
  vardist = vardistExpandParam(vardist, params(kern1.nParams+kern2.nParams+1:kern1.nParams+kern2.nParams+vardist.nParams));
  xx2 = params(kern1.nParams+kern2.nParams+vardist.nParams+1:end);
  x2 = reshape(xx2,[numData/2 numIn]);
  Lminus(i) = full(sum(sum(kernkernVardistPsi2Compute(kern1, kern2, vardist, x2))));
end
params = origParams;
gLDiff = .5*(Lplus - Lminus)/epsilon;
g = [gKern1 gKern2 gVarmeans gVarcovars gInd];
% check firstly the kernel hyperparameters 
kerndiff1 = abs(g(1:kern1.nParams) - gLDiff(1:kern1.nParams));
[g(1:kern1.nParams); gLDiff(1:kern1.nParams)]
pause

index = [kern1.nParams+1:kern1.nParams+kern2.nParams];
kerndiff2 = abs(g(kern1.nParams+1:kern1.nParams+kern2.nParams) ...
            - gLDiff(kern1.nParams+1:kern1.nParams+kern2.nParams));
[g(kern1.nParams+1:kern1.nParams+kern2.nParams); gLDiff(kern1.nParams+1:kern1.nParams+kern2.nParams)]
pause

index = [kern1.nParams+kern2.nParams+1:kern1.nParams+kern2.nParams+(vardist.nParams/2)];
varmeansdiff =  abs(g(index) - gLDiff(index));
[g(index); gLDiff(index)]
pause

index = [kern1.nParams+kern2.nParams+(vardist.nParams/2)+1:kern1.nParams+kern2.nParams+vardist.nParams];
varcovarsdiff =  abs(g(index) - gLDiff(index));
[g(index); gLDiff(index)]
pause
index = [(kern1.nParams+kern2.nParams+vardist.nParams+1):size(g,2)]; 
varinddiff =  abs(g(index) - gLDiff(index)); 
[g(index); gLDiff(index)]
pause
fprintf('Kernel1 hyps max diff: %2.6g.\n', max(kerndiff1));
fprintf('Kernel2 hyps max diff: %2.6g.\n', max(kerndiff2));
fprintf('Variational means max diff: %2.6g.\n', max(varmeansdiff));
fprintf('Variational covars max diff: %2.6g.\n', max(varcovarsdiff));
fprintf('Inducing inputs max diff: %2.6g.\n', max(varinddiff));


