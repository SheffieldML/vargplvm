function [gKern, gVarmeans, gVarcovars] = whiteVardistPsi0Gradient(whitekern, vardist, covGrad)

% WHITEVARDISTPSI0GRADIENT one line description
% FORMAT
% DESC description here.
% RETURN gKern :
% RETURN gVarmeans :
% RETURN gVarcovars :
% RETURN gInd :
% ARG whitekern : the kernel structure associated with the white kernel.
% ARG vardist :
% ARG Z : 
% ARG covGrad : 
%
%
% SEEALSO : others
%
%
% COPYRIGHT : Michalis K. Titsias, 2009
%

% VARGPLVM

  
% the "white" kernel only affects the K_uu matrix (jitter inducing variables)"
% that's why the gKern = 0 
gKern = 0; % covGrad*vardist.numData;
%gKern = covGrad*vardist.numData;
gVarmeans = zeros(1,prod(size(vardist.means))); 
gVarcovars = zeros(1,prod(size(vardist.means))); 