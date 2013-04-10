function [gKern1, gKern2, gVarmeans, gVarcovars, gInd] = rbfard2linard2VardistPsi2Gradient(rbfardKern, linardKern, vardist, Z, covGrad, learnInducing)

% RBFARD2LINARD2VARDISTPSI2GRADIENT description
  
% VARGPLVM

if nargin < 6
    learnInducing = 1;
end

% variational means
N = size(vardist.means,1);
%  inducing variables 
[M Q] = size(Z); 

[Psi2, Pnovar, Psi1] = rbfard2linard2VardistPsi2Compute(rbfardKern, linardKern, vardist, Z);

%-- New: preallocation
    gVarmeans = zeros(N,Q);
    gVarcovars = zeros(N,Q);
    gInd = zeros(M,Q);
%---

% inverse variances
A1 = rbfardKern.inputScales;
A2 = linardKern.inputScales;

% gradient wrt variance of the rbfard2 kernel 
gKernvar = sum(sum(Pnovar.*covGrad));  

AA1 = ones(N,1)*A1;
AS = AA1.*vardist.covars + 1;
SAS = vardist.covars./AS;
MAS = vardist.means./AS;

Zones = ones(size(Z)); 
Kern2tmp = zeros(M,M,Q);  
for n=1:N     
  %
  %
  S_n = vardist.covars(n,:);    
  Mu_n = vardist.means(n,:); 
  AS =  A1.*vardist.covars(n,:) + 1;
  
  ok =  Z*sparse(diag((A1.*A2.*S_n)./AS))*Z' + Zones*sparse(diag((A2.*Mu_n)./AS))*Z'; 
    
  for q=1:vardist.latentDimension 
     %
   
     MZ = (Mu_n(q) - Z(:,q))';
     
     B = (MZ.^2)./AS(q);       
     DerPsi1 = (Psi1(n,:)/AS(q)).*(A1(q)*B - 1);
     DerPsi1 = 0.5*A1(q)*DerPsi1;
       
     tmp = repmat(DerPsi1',[1 M]).*ok ...
                       - repmat(Psi1(n,:),[M 1]).*(Z(:,q)*MZ)*((A1(q)*A2(q))/(AS(q)*AS(q)));   
               
     gVarcovars(n,q) = sum(sum(covGrad.*tmp)); 
     
    
     % kern1 (rbfard2 kernel) hyperprameters 
     B1 = (S_n(q) + B)./AS(q);
     DerPsi1 = -0.5*Psi1(n,:).*B1; 
     Kern2tmp(:,:,q) = Kern2tmp(:,:,q) + repmat(DerPsi1',[1 M]).*ok ...
                         - repmat(Psi1(n,:),[M 1]).*(Z(:,q)*MZ)*((S_n(q)*A2(q))/(AS(q)*AS(q))); 
      
     
     % variational means
     B = (Psi1(n,:).*MZ)./AS(q);
     DerPsi1 = -A1(q)*B; 
    
     tmp = repmat(DerPsi1',[1 M]).*ok + (Psi1(n,:)'*Z(:,q)')*((A2(q))/AS(q));        
     gVarmeans(n,q) = sum(sum(covGrad.*tmp));
     
     %
     %
end
end 


%%%%%%    
for q=1:Q
    %
    AAS = A1(q)*vardist.covars(:,q) + 1;
    S_q = vardist.covars(:,q);
    Mu_q = vardist.means(:,q);
  
    ZZ = Z(:,q)*Z(:,q)';
    
    tmpKern2 = A1(q)*Psi1'*(SAS(:,q)*ones(1,M));
    
    % gradient of kern1 (rbfard2 kernel) hyperprameters
    gKern1(q) = sum(sum(covGrad.*Kern2tmp(:,:,q))); 
     
    % gradient of kern2 (linard2 kernel) hyperprameters
    gKern2(q) = sum(sum( (Psi1'*(MAS(:,q)*Z(:,q)') + tmpKern2.*ZZ).*covGrad ));
    
    
    MZ = repmat(vardist.means(:,q),1,M) - repmat(Z(:,q)',N,1);
    B = Psi1.*(MZ./repmat(AAS,1,M));
    DerPsi1 = A1(q)*B;   

    Q1term = ( Psi1'*(S_q./AAS) )*Z(:,q)'*A1(q)*A2(q);  
    Q2term = sum( ( Psi1./repmat(AAS,1,M) ).*( repmat(Z(:,q)',N,1).*repmat(S_q,1,M)*A1(q)*A2(q) + repmat(Mu_q,1,M)*A2(q) ), 1);
    %
    if learnInducing
        for m=1:M
            %
            Bterm = (vardist.covars./(AA1.*vardist.covars + 1))*sparse(diag(Z(m,:).*A1.*A2))*Z'...
                + (vardist.means./(AA1.*vardist.covars + 1))*sparse(diag(A2))*Z';
            
            Indtmp =  sum(repmat(DerPsi1(:,m),1,M).*Bterm ,1) + Q1term(m,:) + Q2term;
            
            gInd(m,q) = sum(sum(covGrad(m,:).*Indtmp));
            %
        end
    end
    %
end
        

gKern1 = 2*[gKernvar gKern1(:)'];
gKern2 = 2*gKern2(:)';

% gVarmeans is N x Q matrix (N:number of data, Q:latent dimension)
% this will unfold this matrix column-wise
gVarmeans = 2*gVarmeans(:)';

% gVarcovars is N x Q matrix (N:number of data, Q:latent dimension)
% this will unfold this matrix column-wise 
gVarcovars = 2*gVarcovars(:)';

% gInd is M x Q matrix (M:number of inducing variables, Q:latent dimension)
% this will unfold this matrix column-wise 
gInd = 2*gInd(:)'; 


