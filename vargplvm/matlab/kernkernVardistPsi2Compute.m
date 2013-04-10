function [Psi2, Pnovar, Psi1] = kernkernVardistPsi2Compute(kern1, kern2, vardist, Z)

% KERNKERNVARDISTPSI2COMPUTE description.

% VARGPLVM

%

if strcmp(kern2.type, kern1.type)
    error('Kern1 and Kern2 cannot be the same')
    return; 
elseif ~strcmp(kern1.type, 'white') & ~strcmp(kern2.type, 'white')
    %
    swap = 0;
    if strcmp(kern2.type, 'rbfard2') 
       swap = 1;
    elseif strcmp(kern2.type, 'linard2') & ~strcmp(kern1.type, 'rbfard2') 
       swap = 1;   
    end
    %
    if swap == 0
       fhandle = str2func([kern1.type kern2.type 'VardistPsi2Compute']);
       [Psi2, Pnovar, Psi1] = fhandle(kern1, kern2, vardist, Z);
    else
       fhandle = str2func([kern2.type kern1.type 'VardistPsi2Compute']);
       [Psi2, Pnovar, Psi1] = fhandle(kern2, kern1, vardist, Z);
    end    
   %
else 
   %
   % the white kernel gives zeros
   Psi2 = zeros(size(Z,1),size(Z,1)); 
   Pnovar = zeros(size(vardist.means));
   Psi1 = zeros(size(vardist.means));
   %
end


