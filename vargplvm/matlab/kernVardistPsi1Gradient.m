function [gKern, gVarmeans, gVarcovars, gInd] = kernVardistPsi1Gradient(kern, vardist, Z, covGrad, learnInducing)

% KERNVARDISTPSI1GRADIENT description.  

% VARGPLVM

if nargin < 5
    learnInducing = 1;
end

    if ~strcmp(kern.type,'cmpnd')
       % 
       fhandle = str2func([kern.type 'VardistPsi1Gradient']);
       [gKern, gVarmeans, gVarcovars, gInd] = fhandle(kern, vardist, Z, covGrad, learnInducing);    

       % Transformations
       gKern = paramTransformPsi1(kern, gKern); 
       %
    else % the kernel is cmpnd
       %
       fhandle = str2func([kern.comp{1}.type 'VardistPsi1Gradient']);
       [gKern, gVarmeans, gVarcovars, gInd] = fhandle(kern.comp{1}, vardist, Z, covGrad, learnInducing);
       % Transformations
       gKern = paramTransformPsi1(kern.comp{1}, gKern); 
       %
       for i = 2:length(kern.comp)
           %
           fhandle = str2func([kern.comp{i}.type 'VardistPsi1Gradient']);
           [gKerni, gVarmeansi, gVarcovarsi, gIndi] = fhandle(kern.comp{i}, vardist, Z, covGrad, learnInducing);

           % Transformations
           gKerni = paramTransformPsi1(kern.comp{i}, gKerni); 

           gVarmeans = gVarmeans + gVarmeansi; 
           gVarcovars = gVarcovars + gVarcovarsi;
           gInd = gInd + gIndi;
           gKern = [gKern gKerni]; 
           %
       end
       % 
    end

end

% variational variances are positive (This should rather go to
% vargplvmLogLikeGradients)
%gVarcovars = (gVarcovars(:).*vardist.covars(:))'; 



%-----------------------------------------------------
% This applies transformations 
% This must be done similar to kernGradient at some point 
% function gKern = paramTransformPsi1(kern, gKern)
% %
% % 
% if strcmp(kern.type,'rbfard2') | strcmp(kern.type,'rbfardjit')
%    gKern(1) = gKern(1)*kern.variance;
%    gKern(2:end) = gKern(2:end).*kern.inputScales; 
% elseif strcmp(kern.type,'linard2')
%    gKern(1:end) = gKern(1:end).*kern.inputScales;
% elseif strcmp(kern.type,'bias')
%    gKern = gKern*kern.variance;   
% else
%    % do nothing
% end

function gKern = paramTransformPsi1(kern, gKern)

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



