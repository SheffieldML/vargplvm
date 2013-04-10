function [K, P] = linard2VardistPsi1Compute(linard2kern, vardist, Z)

% LINARD2VARDISTPSI1COMPUTE description.

% VARGPLVM


K = kernCompute(linard2kern,vardist.means,Z);

P = [];