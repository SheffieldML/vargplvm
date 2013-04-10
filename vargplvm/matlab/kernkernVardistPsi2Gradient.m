function [gKern1, gKern2, gVarmeans, gVarcovars, gInd] = kernkernVardistPsi2Gradient(kern1, kern2, vardist, Z, covGrad, learnInducing)

% KERNKERNVARDISTPSI2GRADIENT Description
  
% VARGPLVM

if nargin < 6
    learnInducing = true;
end

if strcmp(kern2.type, kern1.type)
    error('Kern1 and Kern2 cannot be the same')
    return; 
elseif ~strcmp(kern1.type, 'white') &  ~strcmp(kern2.type, 'white')
    %
    swap = 0;
    if strcmp(kern2.type, 'rbfard2') 
       swap = 1;
    elseif strcmp(kern2.type, 'linard2') & ~strcmp(kern1.type, 'rbfard2') 
       swap = 1;   
    end
    %
    if swap == 0
       fhandle = str2func([kern1.type kern2.type 'VardistPsi2Gradient']);
       [gKern1, gKern2, gVarmeans, gVarcovars, gInd] = fhandle(kern1, kern2, vardist, Z, covGrad, learnInducing);
    else
       fhandle = str2func([kern2.type kern1.type 'VardistPsi2Gradient']);
       [gKern2, gKern1, gVarmeans, gVarcovars, gInd] = fhandle(kern2, kern1, vardist, Z, covGrad, learnInducing);
    end
    % !!! Do not take transformation here: they applied in the kernkernVardistPsi2Gradient
    %gKern1 = paramTransformPsi2(kern1, gKern1);
    %gKern2 = paramTransformPsi2(kern2, gKern2);
    % variational variances are positive  
    %gVarcovars = (gVarcovars(:).*vardist.covars(:))';
   %
else
   %
    % the white kernel gives zeros
    gKern1 = zeros(1,kern1.nParams); 
    gKern2 = zeros(1,kern2.nParams); 
    gVarmeans = zeros(1,prod(size(vardist.means))); 
    gInd = zeros(1,prod(size(Z))); 
    gVarcovars = zeros(1,prod(size(vardist.covars)));
   %
end


%%
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
