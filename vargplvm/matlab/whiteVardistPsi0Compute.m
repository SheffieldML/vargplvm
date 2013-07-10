function Psi0 = whiteVardistPsi0Compute(whitekern, vardist)

% WHITEVARDISTPSI0COMPUTE one line description
% FORMAT
% DESC description
% RETURN Psi0 : description
% ARG whiteKern : the kernel structure associated with the white kernel.
% ARG vardist : description
%
%
% SEEALSO : others
%
%
% COPYRIGHT : Michalis K. Titsias, 2009
%

% VARGPLVM

% the "white" kernel only affects the K_uu matrix (jitter inducing variables)"
% that's why the psi0_white = 0 
Psi0 = 0; % vardist.numData*whitekern.variance; 
%Psi0 =  vardist.numData*whitekern.variance;
