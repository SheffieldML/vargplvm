function [K, outKern, sumKern, Kgvar] = rbfard2VardistPsi2Compute(rbfardKern, vardist, Z)

% RBFARD2VARDISTPSI2COMPUTE one line description
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

try
    pool_open = matlabpool('size')>0;
catch e
    pool_open = 0;
end

% The conditions for the parallel code to run, is the workers pool to be
% open, the parallel flag to be active and the number of datapoints N to be
% larger than a reasonable threshold (otherwise there is unecessary
% thread-communication overhead).
if pool_open && (isfield(vardist,'parallel') && vardist.parallel) && size(vardist.means,1) > 15
    [K, outKern, sumKern, Kgvar] = rbfard2VardistPsi2ComputePar(rbfardKern, vardist, Z);
else
    [K, outKern, sumKern, Kgvar] = rbfard2VardistPsi2ComputeOrig(rbfardKern, vardist, Z);
end




function [K, outKern, sumKern, Kgvar] = rbfard2VardistPsi2ComputePar(rbfardKern, vardist, Z)

% variational means
N  = size(vardist.means,1);
%  inducing variables 
M = size(Z,1); 

A = rbfardKern.inputScales;
    
% first way
sumKern = zeros(M,M); 

parfor n=1:N
    %    
    AS_n = (1 + 2*A.*vardist.covars(n,:)).^0.5;  
    
    normfactor =  1./prod(AS_n);
    
    %Z_n = (repmat(vardist.means(n,:),[M 1]) - Z)*0.5; 
    Z_n = bsxfun(@minus, vardist.means(n,:), Z)*0.5;
    %Z_n = Z_n.*repmat(sqrt(A)./AS_n,[M 1]);
    Z_n = bsxfun(@times, Z_n, sqrt(A)./AS_n);
    distZ = dist2(Z_n,-Z_n); 
    
    sumKern = sumKern + normfactor*exp(-distZ);  
    %
end
    
% ZZ = Z.*(repmat(sqrt(A),[M 1]));
ZZ =  bsxfun(@times, Z, sqrt(A));
distZZ = dist2(ZZ,ZZ);
outKern = exp(-0.25*distZZ);

Kgvar = rbfardKern.variance*(outKern.*sumKern); 
K = rbfardKern.variance*Kgvar;


function [K, outKern, sumKern, Kgvar] = rbfard2VardistPsi2ComputeOrig(rbfardKern, vardist, Z)

% variational means
N  = size(vardist.means,1);
%  inducing variables 
M = size(Z,1); 

A = rbfardKern.inputScales;
    
% first way
sumKern = zeros(M,M); 
for n=1:N
    %    
    AS_n = (1 + 2*A.*vardist.covars(n,:)).^0.5;  
    
    normfactor =  1./prod(AS_n);
    
    %Z_n = (repmat(vardist.means(n,:),[M 1]) - Z)*0.5; 
    Z_n = bsxfun(@minus, vardist.means(n,:), Z)*0.5;
    %Z_n = Z_n.*repmat(sqrt(A)./AS_n,[M 1]);
    Z_n = bsxfun(@times, Z_n, sqrt(A)./AS_n);
    distZ = dist2(Z_n,-Z_n); 
    
    sumKern = sumKern + normfactor*exp(-distZ);  
    %
end
    
% ZZ = Z.*(repmat(sqrt(A),[M 1]));
ZZ =  bsxfun(@times, Z, sqrt(A));
distZZ = dist2(ZZ,ZZ);
outKern = exp(-0.25*distZZ);

Kgvar = rbfardKern.variance*(outKern.*sumKern); 
K = rbfardKern.variance*Kgvar;

%{
% second way
%sumKern2 = zeros(M,M); 
%for n=1:N
%    norm = [];
%    sumK = zeros(M,M); 
%    for q=1:size(Z,2)
%       norm(q) = sqrt(1 + 2*A(q).*vardist.covars(n,q));
%       
%       ZZ_q = 0.5*(repmat(Z(:,q),[1 M]) + repmat(Z(:,q)',[M 1]));
%       sumK = sumK + (A(q)/(1 + 2*A(q)*vardist.covars(n,q)))*((vardist.means(n,q) - ZZ_q).^2);
%    end
%    sumKern2 = sumKern2 + exp(-sumK)/prod(norm); 
%end
%(1/prod(norm))
%
%outKern2 = zeros(M,M); 
%for q=1:size(Z,2)
%    outKern2 = outKern2 + A(q)*dist2(Z(:,q),Z(:,q)); 
%end
%outKern2 = exp(-0.25*outKern2); 
%sum(sum(abs(sumKern - sumKern2)))
%sum(sum(abs(outKern - outKern2)))
%Kgvar = rbfardKern.variance*(outKern2.*sumKern2); 
%K = rbfardKern.variance*Kgvar;

% third way 
%zmAzm = sum((Z.*Z).*repmat(A,[M 1]),2);
%Azm = Z.*repmat(A,[M 1]);
%outKern2 = repmat(zmAzm,[1 M]) + repmat(zmAzm',[M 1]);
%outKern2 = (rbfardKern.variance^2)*exp(-0.5*outKern2); 
%
%sumKern2 = zeros(M,M); 
%for n=1:N
%    %    
%    AS_n = (1 + 2*A.*vardist.covars(n,:)).^0.5;  
%    
%    mu_n = vardist.means(n,:);
%    mmun = mu_n*(inv(diag(vardist.covars(n,:)))*mu_n'); 
%    
%    normfactor =  1./prod(AS_n);
%    for m=1:M
%    for mp=1:M 
%        argExp = Z(m,:).*A + Z(mp,:).*A + vardist.means(n,:)./vardist.covars(n,:); 
%        argExp = argExp*(inv(2*diag(A) + diag(1./vardist.covars(n,:)))*argExp'); 
%        tmpKern(m,mp) = exp(-0.5*(mmun) +0.5*argExp); 
%    end
%    end
%    sumKern2 = sumKern2 + normfactor*tmpKern;
%    %
%end
%
%KK = outKern2.*sumKern2;
%sum(sum(abs(KK - K)))
%K
%KK
%pause
%sum(sum(abs(outKern - outKern2)))
 %}   
