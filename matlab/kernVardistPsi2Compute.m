function [Psi2 P] = kernVardistPsi2Compute(kern, vardist, Z)

% KERNVARDISTPSI2COMPUTE description.  

% VARGPLVM

% compute first the "square" terms of Psi2 
if ~strcmp(kern.type,'cmpnd')
  %   
  fhandle = str2func([kern.type 'VardistPsi2Compute']);
  [Psi2 P] = fhandle(kern, vardist, Z);
  %  
else % the kernel is cmpnd
  %
  fhandle = str2func([kern.comp{1}.type 'VardistPsi2Compute']);
  [Psi2 P] = fhandle(kern.comp{1}, vardist, Z);
  %
  for i = 2:length(kern.comp)
      %
      fhandle = str2func([kern.comp{i}.type 'VardistPsi2Compute']);
      [Ptmp P] = fhandle(kern.comp{i}, vardist, Z);
      Psi2 = Psi2 + Ptmp;
      %
  end
  %
end

% compute the cross-kernel terms of Psi2
if strcmp(kern.type,'cmpnd') & (length(kern.comp)>1)
  for i=1:length(kern.comp)
    for j=i+1:length(kern.comp)
       %[i j]
       [Ptmp, P] = kernkernVardistPsi2Compute(kern.comp{i}, kern.comp{j}, vardist, Z); 
       Psi2 = Psi2 + Ptmp;
   end
  end
end

