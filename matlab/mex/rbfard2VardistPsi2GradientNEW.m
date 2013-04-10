function [gKern, gVarmeans, gVarcovars, gInd] = rbfard2VardistPsi2GradientNEW(rbfardKern, vardist, Z, covGrad)

% RBFARD2VARDISTPSI2GRADIENTNEW description.
  
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
    [gKern, gVarmeans, gVarcovars, gInd] = rbfard2VardistPsi2GradientPar(rbfardKern, vardist, Z, covGrad);
else
    [gKern, gVarmeans, gVarcovars, gInd] = rbfard2VardistPsi2GradientOrig(rbfardKern, vardist, Z, covGrad);
end




function [gKern, gVarmeans, gVarcovars, gInd] = rbfard2VardistPsi2GradientPar(rbfardKern, vardist, Z, covGrad)
% variational means
N = size(vardist.means,1);
%  inducing variables 
[M Q] = size(Z); 


% evaluate the kernel matrix 
[K, outKern, sumKern, Kgvar] = rbfard2VardistPsi2Compute(rbfardKern, vardist, Z);

% inverse variances
A = rbfardKern.inputScales;

% gradient wrt variance of the kernel 
gKernvar = 2*sum(sum(Kgvar.*covGrad));  


% 1) line compute 0.5*(z_mq + z_m'q) for any q and store the result in a "M x Q x M" 
%  matrix where M is the number of inducing points and Q the latent dimension
% 2) line compute the z_mq - z_m'q, for any q
ZmZm  = zeros(M,Q,M);
ZmDZm = zeros(M,Q,M);
for q=1:size(Z,2)
  ZmZm(:,q,:) = 0.5*(repmat(Z(:,q),[1 1 M]) + repmat(reshape(Z(:,q),[1 1 M]),[M 1 1]));
  ZmDZm(:,q,:) = repmat(Z(:,q),[1 1 M]) - repmat(reshape(Z(:,q),[1 1 M]),[M 1 1]);
end

% compute the terms 2 a_q s_nq^2 + 1, for n and q and srore the result in a 
% "N x Q" matrix
asPlus1 = 2*(repmat(A,[N 1]).*vardist.covars) + 1; 
% compute the terms a_q/(2 a_q s_nq^2 + 1), for n and q and store the result in a 
% "N x Q" matrix
aDasPlus1 = repmat(A,[N 1])./asPlus1; 

covGrad = (rbfardKern.variance^2)*(covGrad.*outKern);
covGrad = reshape(covGrad,[M 1 M]);
sumKern = reshape(sumKern,[M 1 M]);
Amq = repmat(A,[M 1]);
partInd1 = - Amq.*sum(ZmDZm.*repmat(sumKern.*covGrad,[1 Q 1]),3);
partInd2 = zeros(M,Q);

partA1 = - 0.25*sum(sum((ZmDZm.*ZmDZm).*repmat(sumKern.*covGrad,[1 Q 1]),3),1);
partA2 = zeros(1,Q);

gVarcovars = zeros(N,Q); 
gVarmeans = zeros(N,Q);

% Compute the gradient wrt lengthscales, variational means and variational variances  
% For loop over training points  
parfor n=1:N
    %
    %  
    mu_n = vardist.means(n,:); 
    s2_n = vardist.covars(n,:); 
    AS_n = asPlus1(n,:);  
     
    %MunZmZm = repmat(mu_n, [M 1 M]) - ZmZm;
    MunZmZm = bsxfun(@minus,mu_n,ZmZm);
    %MunZmZmA = MunZmZm./repmat(AS_n,[M 1 M]);
    MunZmZmA =  bsxfun(@rdivide, MunZmZm, AS_n);
    
    %k2Kern_n = sum((MunZmZm.^2).*repmat(aDasPlus1(n,:),[M 1 M]),2);    
    k2Kern_n = sum(  bsxfun(@times, MunZmZm.^2,aDasPlus1(n,:)),2);
    
    k2Kern_n = exp(-k2Kern_n)/prod(sqrt(AS_n));
    
    % derivatives wrt to variational means
    k2ncovG = repmat(k2Kern_n.*covGrad,[1 Q 1]); 
    %tmp2 = tmp + reshape(diag(diag(squeeze(tmp))),[M 1 M]);
    %diagCorr = diag(diag(squeeze(tmp))); 
    tmp = MunZmZmA.*k2ncovG;
    tmp = sum(tmp,3);
    gVarmeans(n,:) = - 2*A.*(sum(tmp,1));
    
    % derivatives wrt inducing inputs 
    %diagCorr = diagCorr*(repmat(mu_n,[M 1]) - Z).*repmat(aDasPlus1(n,:),[M 1]);
    %partInd2 = partInd2 + Amq.*(sum(tmp,3) + diagCorr);
    partInd2 = partInd2 + Amq.*tmp;
    
    
    % Derivative wrt input scales  
    MunZmZmA = MunZmZmA.*MunZmZm; 
    %partA2 = partA2 + sum(sum(((MunZmZmA + repmat(s2_n,[M 1 M])).*k2ncovG)./repmat(AS_n,[M 1 M]),1),3);
    tmppartA2 = bsxfun(@plus, MunZmZmA,s2_n).*k2ncovG;
    partA2 = partA2 + sum(sum( bsxfun(@rdivide, tmppartA2, AS_n), 1),3);
    
    % derivatives wrt variational diagonal covariances 
    %MunZmZmA = MunZmZmA.*repmat(A,[M 1 M]);
    MunZmZmA = bsxfun(@times, MunZmZmA, A);
    %gVarcovars(n,:) = sum(sum(repmat(aDasPlus1(n,:),[M 1 M]).*(2*MunZmZmA - 1).*k2ncovG,1),3);
    gVarcovars(n,:) = sum(sum( bsxfun(@times, (2*MunZmZmA - 1).*k2ncovG, aDasPlus1(n,:)),1),3);
    
    %ZmZm1 = k2kernCompute(A, mu_n, cov_n, Z); 
    %
    %AS_n = (1 + 2*A.*vardist.covars(n,:)).^0.5;  
    %
    %normfactor =  1./prod(AS_n);
    %
    %Z_n = (repmat(vardist.means(n,:),[M 1]) - Z)*0.5; 
    %Z_n = Z_n.*repmat(sqrt(A)./AS_n,[M 1]);
    %distZ = dist2(Z_n,-Z_n); 
    %
    %sumKern = sumKern + normfactor*exp(-distZ);  
    %
end

gInd = partInd1 + 2*partInd2; 

gKernlengcs = partA1 - partA2; 
gKern = [gKernvar gKernlengcs];

% gVarmeans is N x Q matrix (N:number of data, Q:latent dimension)
% this will unfold this matrix column-wise 
%gVarmeans = gVarmeans'; 
gVarmeans = gVarmeans(:)'; 

% gVarcovars is N x Q matrix (N:number of data, Q:latent dimension)
% this will unfold this matrix column-wise 
%gVarcovars = gVarcovars'; 
gVarcovars = gVarcovars(:)';

% gInd is M x Q matrix (M:number of inducing variables, Q:latent dimension)
% this will unfold this matrix column-wise 
%gInd = gInd'; 
gInd = gInd(:)'; 





function [gKern, gVarmeans, gVarcovars, gInd] = rbfard2VardistPsi2GradientOrig(rbfardKern, vardist, Z, covGrad)
% variational means
N = size(vardist.means,1);
%  inducing variables 
[M Q] = size(Z); 


% evaluate the kernel matrix 
[K, outKern, sumKern, Kgvar] = rbfard2VardistPsi2Compute(rbfardKern, vardist, Z);

% inverse variances
A = rbfardKern.inputScales;

% gradient wrt variance of the kernel 
gKernvar = 2*sum(sum(Kgvar.*covGrad));  


% 1) line compute 0.5*(z_mq + z_m'q) for any q and store the result in a "M x Q x M" 
%  matrix where M is the number of inducing points and Q the latent dimension
% 2) line compute the z_mq - z_m'q, for any q
ZmZm  = zeros(M,Q,M);
ZmDZm = zeros(M,Q,M);
for q=1:size(Z,2)
  ZmZm(:,q,:) = 0.5*(repmat(Z(:,q),[1 1 M]) + repmat(reshape(Z(:,q),[1 1 M]),[M 1 1]));
  ZmDZm(:,q,:) = repmat(Z(:,q),[1 1 M]) - repmat(reshape(Z(:,q),[1 1 M]),[M 1 1]);
end

% compute the terms 2 a_q s_nq^2 + 1, for n and q and srore the result in a 
% "N x Q" matrix
asPlus1 = 2*(repmat(A,[N 1]).*vardist.covars) + 1; 
% compute the terms a_q/(2 a_q s_nq^2 + 1), for n and q and store the result in a 
% "N x Q" matrix
aDasPlus1 = repmat(A,[N 1])./asPlus1; 

covGrad = (rbfardKern.variance^2)*(covGrad.*outKern);
covGrad = reshape(covGrad,[M 1 M]);
sumKern = reshape(sumKern,[M 1 M]);
Amq = repmat(A,[M 1]);
partInd1 = - Amq.*sum(ZmDZm.*repmat(sumKern.*covGrad,[1 Q 1]),3);
partInd2 = zeros(M,Q);

partA1 = - 0.25*sum(sum((ZmDZm.*ZmDZm).*repmat(sumKern.*covGrad,[1 Q 1]),3),1);
partA2 = zeros(1,Q);

gVarcovars = zeros(N,Q); 
gVarmeans = zeros(N,Q);

% Compute the gradient wrt lengthscales, variational means and variational variances  
% For loop over training points  
% for n=1:N
%     %
%     %  
%     mu_n = vardist.means(n,:); 
%     s2_n = vardist.covars(n,:); 
%     AS_n = asPlus1(n,:);  
%      
%     %MunZmZm = repmat(mu_n, [M 1 M]) - ZmZm;
%     MunZmZm = bsxfun(@minus,mu_n,ZmZm);
%     %MunZmZmA = MunZmZm./repmat(AS_n,[M 1 M]);
%     MunZmZmA =  bsxfun(@rdivide, MunZmZm, AS_n);
%     
%     %k2Kern_n = sum((MunZmZm.^2).*repmat(aDasPlus1(n,:),[M 1 M]),2);    
%     k2Kern_n = sum(  bsxfun(@times, MunZmZm.^2,aDasPlus1(n,:)),2);
%     
%     k2Kern_n = exp(-k2Kern_n)/prod(sqrt(AS_n));
%     
%     % derivatives wrt to variational means
%     k2ncovG = repmat(k2Kern_n.*covGrad,[1 Q 1]); 
%     %tmp2 = tmp + reshape(diag(diag(squeeze(tmp))),[M 1 M]);
%     %diagCorr = diag(diag(squeeze(tmp))); 
%     tmp = MunZmZmA.*k2ncovG;
%     tmp = sum(tmp,3);
%     gVarmeans(n,:) = - 2*A.*(sum(tmp,1));
%     
%     % derivatives wrt inducing inputs 
%     %diagCorr = diagCorr*(repmat(mu_n,[M 1]) - Z).*repmat(aDasPlus1(n,:),[M 1]);
%     %partInd2 = partInd2 + Amq.*(sum(tmp,3) + diagCorr);
%     partInd2 = partInd2 + Amq.*tmp;
%     
%     
%     % Derivative wrt input scales  
%     MunZmZmA = MunZmZmA.*MunZmZm; 
%     %partA2 = partA2 + sum(sum(((MunZmZmA + repmat(s2_n,[M 1 M])).*k2ncovG)./repmat(AS_n,[M 1 M]),1),3);
%     tmppartA2 = bsxfun(@plus, MunZmZmA,s2_n).*k2ncovG;
%     partA2 = partA2 + sum(sum( bsxfun(@rdivide, tmppartA2, AS_n), 1),3);
%     
%     % derivatives wrt variational diagonal covariances 
%     %MunZmZmA = MunZmZmA.*repmat(A,[M 1 M]);
%     MunZmZmA = bsxfun(@times, MunZmZmA, A);
%     %gVarcovars(n,:) = sum(sum(repmat(aDasPlus1(n,:),[M 1 M]).*(2*MunZmZmA - 1).*k2ncovG,1),3);
%     gVarcovars(n,:) = sum(sum( bsxfun(@times, (2*MunZmZmA - 1).*k2ncovG, aDasPlus1(n,:)),1),3);
%     
%     %ZmZm1 = k2kernCompute(A, mu_n, cov_n, Z); 
%     %
%     %AS_n = (1 + 2*A.*vardist.covars(n,:)).^0.5;  
%     %
%     %normfactor =  1./prod(AS_n);
%     %
%     %Z_n = (repmat(vardist.means(n,:),[M 1]) - Z)*0.5; 
%     %Z_n = Z_n.*repmat(sqrt(A)./AS_n,[M 1]);
%     %distZ = dist2(Z_n,-Z_n); 
%     %
%     %sumKern = sumKern + normfactor*exp(-distZ);  
%     %
% end

%%%%%%%%5
%writeMatrixFile('means.txt', vardist.means);
%writeMatrixFile('covars.txt',vardist.covars);
%writeMatrixFile('asPlus1.txt',asPlus1);
%writeMatrixFile('aDasPlus1.txt', aDasPlus1);
%writeMatrixFile('ZmZm.txt', ZmZm);
%writeMatrixFile('covGrad.txt', covGrad);

covGrad_dim2 = size(covGrad,2);
% mex vargplvm.cpp


%vargplvm([M,N,Q],A,ZmZm_arr, covGrad_dim2)
ZmZmRes = reshape(ZmZm, size(ZmZm,1),size(ZmZm,2)*size(ZmZm,3));
covGradRes = reshape(covGrad, size(covGrad,1), size(covGrad,2)*size(covGrad,3));

[partInd2,partA2,gVarmeans,gVarcovars] = vargplvm([M,N,Q],A,covGrad_dim2, vardist.means', vardist.covars', asPlus1', aDasPlus1', ZmZmRes, covGradRes);

partInd2 = reshape(partInd2, M, Q);
gVarmeans = reshape(gVarmeans, N,Q);
gVarcovars = reshape(gVarcovars, N,Q);


%partInd2Old= dlmread('partInd2.txt',  '\t');
%partInd2Old = partInd2Old';
%partA2Old= dlmread('partA2.txt',  '\t');
%partA2Old = partA2Old';
%gVarmeansOld= dlmread('gVarmeans.txt',  '\t');
%gVarmeansOld = gVarmeansOld';
%gVarcovarsOld= dlmread('gVarcovars.txt', '\t');
%gVarcovarsOld = gVarcovarsOld';
%sum(sum(abs(gVarmeans - gVarmeansOld)))
%sum(sum(abs(gVarcovars - gVarcovarsOld)))
%sum(sum(abs(partInd2 - partInd2Old)))
%%%%%%%


gInd = partInd1 + 2*partInd2; 

gKernlengcs = partA1 - partA2; 
gKern = [gKernvar gKernlengcs];

% gVarmeans is N x Q matrix (N:number of data, Q:latent dimension)
% this will unfold this matrix column-wise 
%gVarmeans = gVarmeans'; 
gVarmeans = gVarmeans(:)'; 

% gVarcovars is N x Q matrix (N:number of data, Q:latent dimension)
% this will unfold this matrix column-wise 
%gVarcovars = gVarcovars'; 
gVarcovars = gVarcovars(:)';

% gInd is M x Q matrix (M:number of inducing variables, Q:latent dimension)
% this will unfold this matrix column-wise 
%gInd = gInd'; 
gInd = gInd(:)'; 




