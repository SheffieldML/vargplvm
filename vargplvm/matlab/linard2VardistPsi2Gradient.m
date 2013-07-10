function [gKern, gVarmeans, gVarcovars, gInd] = linard2VardistPsi2Gradient(linard2Kern, vardist, Z, covGrad, learnInducing)

% LINARD2VARDISTPSI2GRADIENT description.

% VARGPLVM

if nargin < 5
    learnInducing = 1;
end

% inverse variances
A = linard2Kern.inputScales;

Psi1 = linard2VardistPsi1Compute(linard2Kern, vardist, Z);

Amat = sparse(diag(A));
K1 = (vardist.means'*vardist.means + diag(sum(vardist.covars,1)));
K = Amat*K1*Amat;

% TYPICAL WAY
%sumS = sum(vardist.covars,1);
%for q=1:vardist.latentDimension
%   sc = sum(vardist.means(:,q)*ones(1,size(Z,1)).*Psi1,1); 
%   
%   ZZ = Z(:,q)*Z(:,q)';
%      
%   gKern(q)  = sum(sum( (Z(:,q)*sc + (sumS(q)*A(q))*ZZ).*covGrad )); 
%   
%   gVarmeans(:,q) =  A(q)*sum(Psi1*((ones(size(Z,1),1)*Z(:,q)').*covGrad),2);
%   
%   gVarcovars(q) = sum(sum((ZZ.*covGrad)));   
%end
%gKern = 2*gKern(:)';
%gVarmeans = 2*gVarmeans(:)'; 
%gVarcovars = (A.^2).*gVarcovars;
%gVarcovars = ones(size(Psi1,1),1)*gVarcovars;
%gVarcovars = gVarcovars(:)';
%gInd = covGrad*Z*K;
%gInd = 2*gInd(:)';
%%% end of typical way 

% FAST WAY
AZ = Z*Amat;
AZK = AZ*K1;
gKern = sum((covGrad*Z).*AZK,1);
gVarmeans = (Psi1*(covGrad*Z))*Amat;
AZ = AZ*Amat;
gVarcovars = ones(size(Psi1,1),1)*sum((covGrad*Z).*AZ,1);
if learnInducing
    gInd = covGrad*Z*K;
else
    [M Q] = size(Z); 
    gInd = zeros(M, Q);
end
%
gKern = 2*gKern(:)';
gVarmeans = 2*gVarmeans(:)'; 
gVarcovars = gVarcovars(:)';
gInd = 2*gInd(:)';
%%%  end of FAST WAY

%sum(sum(abs(gKern1-gKern)))
%sum(sum(abs(gVarmeans1 - gVarmeans)))
%sum(sum(abs(gVarcovars1 - gVarcovars)))
%sum(sum(abs(gInd1 - gInd)))
%pause
