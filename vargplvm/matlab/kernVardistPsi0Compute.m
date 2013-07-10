function Psi0 = kernVardistPsi0Compute(kern, vardist)

% KERNVARDISTPSI0COMPUTE description.  

% VARGPLVM


if ~strcmp(kern.type,'cmpnd')
  %  
  fhandle = str2func([kern.type 'VardistPsi0Compute']);
  Psi0 = fhandle(kern, vardist);
  %
else  % the kernel is cmpnd
  % 
  fhandle = str2func([kern.comp{1}.type 'VardistPsi0Compute']);
  Psi0 = fhandle(kern.comp{1}, vardist);
  %
  for i = 2:length(kern.comp)
      %
      fhandle = str2func([kern.comp{i}.type 'VardistPsi0Compute']);
      Psi0 = Psi0 + fhandle(kern.comp{i}, vardist);
      %
  end
  %
end

