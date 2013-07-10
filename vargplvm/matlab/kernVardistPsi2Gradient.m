function [gKern, gVarmeans, gVarcovars, gInd] = kernVardistPsi2Gradient(kern, vardist, Z, covGrad, learnInducing)

% KERNVARDISTPSI2GRADIENT description.  

% VARGPLVM

if nargin < 5
    % If learnInducing == 0, then gInd will be returned as a a vector of
    % zeros. It would be better if nargout was checked and gInd was not
    % returned at all if not needed, but the code is very hierarchical and
    % recursive and that would require more complicated code than just
    % including the "learnInducing" flag.
    learnInducing = 1;
end


    %% compute first the "square" terms of Psi2 
    %
    if ~strcmp(kern.type,'cmpnd')
       % 
       fhandle = str2func([kern.type 'VardistPsi2Gradient']);
       [gKern, gVarmeans, gVarcovars, gInd] = fhandle(kern, vardist, Z, covGrad, learnInducing);    

       % Transformations
       gKern = paramTransformPsi2(kern, gKern); 
       %
    else % the kernel is cmpnd
       %
       fhandle = str2func([kern.comp{1}.type 'VardistPsi2Gradient']);
       [gKern, gVarmeans, gVarcovars, gInd] = fhandle(kern.comp{1}, vardist, Z, covGrad, learnInducing);
       % Transformations
       gKern = paramTransformPsi2(kern.comp{1}, gKern); 
       %
       for i = 2:length(kern.comp)
           %
           fhandle = str2func([kern.comp{i}.type 'VardistPsi2Gradient']);
           [gKerni, gVarmeansi, gVarcovarsi, gIndi] = fhandle(kern.comp{i}, vardist, Z, covGrad, learnInducing);

           % Transformations
           gKerni = paramTransformPsi2(kern.comp{i}, gKerni); 

           gVarmeans = gVarmeans + gVarmeansi; 
           gVarcovars = gVarcovars + gVarcovarsi;
           gInd = gInd + gIndi;
           gKern = [gKern gKerni]; 
           %
       end
       % 
    end


    %% compute the cross-kernel terms of Psi2 
    %
    if strcmp(kern.type,'cmpnd') & (length(kern.comp)>1)
      % 
      index{1}=1:kern.comp{1}.nParams; 
      for i=2:length(kern.comp)
           index{i} = (index{i-1}(end)+1):index{i-1}(end)+kern.comp{i}.nParams;
      end
      %
      for i=1:length(kern.comp)  
        for j=i+1:length(kern.comp)

           [gKerni, gKernj, gVarm, gVarc, gI] = kernkernVardistPsi2Gradient(kern.comp{i}, kern.comp{j}, vardist, Z, covGrad, learnInducing);  

           % Transformations
           gKerni = paramTransformPsi2(kern.comp{i}, gKerni);
           gKernj = paramTransformPsi2(kern.comp{j}, gKernj);

           %index{i}
           %index{j}
           %gKerni
           %gKernj
           %gKern
           %pause

           gKern(index{i}) = gKern(index{i}) + gKerni;
           gKern(index{j}) = gKern(index{j}) + gKernj;  

           gVarmeans = gVarmeans + gVarm; 
           gVarcovars = gVarcovars + gVarc;
           gInd = gInd + gI;
           %
       end
      end
    end
end


% variational variances are positive (This should rather go to
% vargplvmLogLikeGradients) 
% gVarcovars = (gVarcovars(:).*vardist.covars(:))';


%-----------------------------------------------------
% This applies transformations 
% This must be done similar to kernGradient at some point 
% function gKern = paramTransformPsi2(kern, gKern)
% %
% %
% 
% if strcmp(kern.type,'rbfard2') | strcmp(kern.type,'rbfardjit')
%    gKern(1) = gKern(1)*kern.variance;
%    gKern(2:end) = gKern(2:end).*kern.inputScales;
% elseif strcmp(kern.type,'linard2')
%    gKern(1:end) = gKern(1:end).*kern.inputScales;
% elseif strcmp(kern.type,'bias')
%    gKern = gKern*kern.variance; 
% end

function gKern = paramTransformPsi2(kern, gKern)

    fhandle = str2func([kern.type 'KernExtractParam']);
    params  = fhandle(kern);

    if length(kern.transforms) > 1
        for i=1:length(kern.transforms)       
            fhandle = str2func([kern.transforms(i).type 'Transform']);
            if ~isfield(kern.transforms(i), 'transformsettings')
                gKern(kern.transforms(i).index) = gKern(kern.transforms(i).index).*fhandle(params(kern.transforms(i).index), 'gradfact');
            else
                gKern(kern.transforms(i).index) = gKern(kern.transforms(i).index).*fhandle(params(kern.transforms(i).index), 'gradfact', kern.transforms(i).transformsettings);
            end
        end
    else
        fhandle = str2func([kern.transforms.type 'Transform']);
        if ~isfield(kern.transforms, 'transformsettings')
            gKern(kern.transforms.index) = gKern(kern.transforms.index).*fhandle(params(kern.transforms.index), 'gradfact');
        else
            gKern(kern.transforms.index) = gKern(kern.transforms.index).*fhandle(params(kern.transforms.index), 'gradfact', kern.transforms.transformsettings);
        end
    end      
    
end
