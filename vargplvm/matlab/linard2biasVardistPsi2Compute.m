function [Psi2, Pnobias, Psi1] = linard2biasVardistPsi2Compute(linardKern, biasKern, vardist, Z)

% LINARD2BIASVARDISTPSI2COMPUTE one line description
% FORMAT
% DESC description
% RETURN Psi2 : description
% RETURN Pnobias : description
% RETURN Psi1 : description
% ARG linard2Kern : the kernel structure associated with the linard2 kernel.
% ARG biasKern : description
% ARG vardist : description
% ARG Z : description
%
%
% SEEALSO : others
%
%
% COPYRIGHT : Michalis K. Titsias, 2009
%
% COPYRIGHT : Neil D. Lawrence, 2009
%

% VARGPLVM


Psi1 = linard2VardistPsi1Compute(linardKern, vardist, Z); 

sumPsi1 = sum(Psi1,1); 

Psi2 = ones(size(Z,1),1)*sumPsi1; 

Pnobias = Psi2 + Psi2'; 
Psi2 = biasKern.variance*Pnobias; 

% second naive way
%Psi22 = zeros(size(Z,1),size(Z,1));
%for j=1:size(Z,1)
%    for i=1:size(Z,1)
%        Psi22(j,i) = biasKern.variance*(sum(Psi1(:,j)) + sum(Psi1(:,i)));
%    end
%end
%sum(sum(abs(Psi2 - Psi22))) 

