function [K, outKern, sumKern, Kgvar] = rbfardjitVardistPsi2Compute(rbfardKern, vardist, Z)

% RBFARDJITVARDISTPSI2COMPUTE one line description
% FORMAT
% DESC description
% RETURN K : description
% RETURN outKern : description
% RETURN sumKern : description
% RETURN Kgvar : description
% ARG rbfardKern : the kernel structure associated with the rbfard2 kernel.
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

[K, outKern, sumKern, Kgvar] = rbfard2VardistPsi2Compute(rbfardKern, vardist, Z);