function [Psi2, Pnovar, Psi1] = rbfard2linard2VardistPsi2Compute(rbfardKern, linardKern, vardist, Z)

% RBFARD2LINARD2VARDISTPSI2COMPUTE description

% VARGPLVM

% variational means
N  = size(vardist.means,1);
%  inducing variables 
M = size(Z,1); 

[Psi1 Pnovar] = rbfard2VardistPsi1Compute(rbfardKern, vardist, Z); 

A1 = rbfardKern.inputScales;
A2 = linardKern.inputScales;

AA1 = ones(N,1)*A1;
AS = AA1.*vardist.covars + 1;
QsumS = AA1.*(vardist.covars./AS);
QsumM = vardist.means./AS;

AZ = sparse(diag(A2))*Z';
Pnovar = (Z.*(Pnovar'*QsumS)  + Pnovar'*QsumM)*AZ;
Psi2 = rbfardKern.variance*(Pnovar + Pnovar');

% naive way
%Psi22 = zeros(M,M);
%S = vardist.covars;  
%Mu = vardist.means;
%for m1=1:M
%for m2=1:M
%   for n=1:N 
%      ok = (S(n,:).*A1.*A2.*Z(m1,:).*Z(m2,:) + Mu(n,:).*A2.*Z(m2,:))./(1 + A1.*S(n,:)); 
%      Psi22(m1,m2) = Psi22(m1,m2) + Psi1(n,m1)*sum(ok);   
%   end
%end
%end
        
%Psi2 = Psi22 + Psi22';

%Psi22 
%Psi2
%sum(sum(abs(Psi2 - Psi22)))
%
%pause
