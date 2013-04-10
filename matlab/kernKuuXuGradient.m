function gInd = kernKuuXuGradient(kern, Z, covGrad)

% KERNKUUXUGRADIENT Description

% VARGPLVM

if ~strcmp(kern.type,'cmpnd')
  % 
  fhandle = str2func([kern.type 'KuuXuGradient']);
  gInd = fhandle(kern, Z, covGrad);    
  %
else % the kernel is cmpnd
     %
  fhandle = str2func([kern.comp{1}.type 'KuuXuGradient']);
  gInd = fhandle(kern.comp{1}, Z, covGrad);
  %
  for i = 2:length(kern.comp)
    %
    fhandle = str2func([kern.comp{i}.type 'KuuXuGradient']);
    gIndi = fhandle(kern.comp{i}, Z, covGrad);
    
    gInd = gInd + gIndi;
    %
  end
  % 
end



%%
function  gInd = linard2KuuXuGradient(linard2Kern, Z, covGrad)
%
  A = linard2Kern.inputScales;
  A = sparse(diag(A));
  ZA = Z*A;
  %M = size(Z,1);
  %   for j=1:size(Z,1)
  %       for q=1:size(Z,2)
  %           ok = zeros(M,M);
  %           ok(:,j) = A(q,q)*Z(:,q); 
  %           %ok(j,:) = A(q,q)*Z(:,q)'; 
  %           %ok(j,j) = 2*ok(j,j);
  %           gInd(j,q) = sum(sum(ok.*covGrad));
  %       end
  %   end
  %   %gInd(:,q) = sum(repmat(ZA(:,q),1,size(Z,1)).*covGrad,2);
  %gInd = 2*gInd(:)';
  gInd = covGrad*ZA;
  gInd = 2*gInd(:)';
%sum(sum(abs(gInd2(:) - gInd(:))))


%%
function  gInd = biasKuuXuGradient(linard2Kern, Z, covGrad)
%
  gInd = zeros(1,prod(size(Z))); 
  
%%
function  gInd = whiteKuuXuGradient(linard2Kern, Z, covGrad)
%
  gInd = zeros(1,prod(size(Z))); 
  
  


%-----------------------------------------------------
% This applies transformations 
% This must be done similar to kernGradient at some point 
function gKern = paramTransformPsi2(kern, gKern)
%
%

if strcmp(kern.type,'rbfard2')
   gKern(1) = gKern(1)*kern.variance;
   gKern(2:end) = gKern(2:end).*kern.inputScales;
elseif strcmp(kern.type,'linard2')
   gKern(1:end) = gKern(1:end).*kern.inputScales;
elseif strcmp(kern.type,'bias')
   gKern = gKern*kern.variance; 
end
