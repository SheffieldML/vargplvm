function [Psi2, P] = linard2VardistPsi2Compute(linard2kern, vardist, Z)

% LINARD2VARDISTPSI2COMPUTE description.

% VARGPLVM


%Psi1 = linard2KernCompute(linard2kern, vardist.means, Z);
%sqrtAS = sparse(diag(linard2kern.inputScales.*sqrt(sum(vardist.covars,1))));
%Zsc = Z*sqrtAS;
%P = Psi1'*Psi1;
%Psi2 = P + Zsc*Zsc';

ZA = Z*sparse(diag(linard2kern.inputScales));

P = vardist.means'*vardist.means + diag(sum(vardist.covars,1)); 

Psi2 = ZA*P*ZA';





