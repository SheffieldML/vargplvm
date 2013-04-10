function model = vargplvmParamInit(model,Y,X, options)

% VARGPLVMPARAMINIT Initialize the variational GPLVM from the data
% COPYRIGHT: Michalis Titsias 2009-2011
% VARGPLVM
  
  
% if input dimension is equal to the latent dimension, 
% then  initialize the variational mean to the normalzed training data
%if model.d == model.q 
%    X = Y; 
%    model.vardist.means = X;  
%    % inducing points 
%    ind = randperm(model.N);
%    ind = ind(1:model.k);
%    model.X_u = X(ind, :);
%end

if nargin < 4
    options = [];
end

if ~strcmp(model.kern.type,'cmpnd')
   % 
   if strcmp(model.kern.type,'rbfard2')  | strcmp(model.kern.type,'rbfardjit') 
      % 
       model.kern.inputScales = 5./(((max(X)-min(X))).^2);
       %model.kern.variance = max(var(Y)); % !!! better use mean(var(Y)) here...
       if isfield(model, 'mOrig')
            model.kern.variance = var(model.mOrig(:)); %mean(var(Y)); % NEW!!!
       else
            model.kern.variance = var(model.m(:)); %mean(var(Y)); % NEW!!!
       end
      %
   elseif strcmp(model.kern.type,'linard2')
      %
      model.kern.inputScales = 5./(((max(X)-min(X))).^2);
      %   
   end
   %
else
   %
   for i = 1:length(model.kern.comp)
      %
      if strcmp(model.kern.comp{i}.type,'rbfard2') | strcmp(model.kern.type,'rbfardjit') 
      % 
         model.kern.comp{i}.inputScales = 5./(((max(X)-min(X))).^2);
         %model.kern.comp{i}.variance = max(var(Y));
         if isfield(model, 'mOrig')
             model.kern.comp{i}.variance = var(model.mOrig(:));
         else
            model.kern.comp{i}.variance = var(model.m(:));
         end
      %
      elseif strcmp(model.kern.comp{i}.type,'linard2')
      %
       model.kern.comp{i}.inputScales = 0.01*max(var(Y))*ones(1,size(X,2));% %5./(((max(X)-min(X))).^2);
      %   
      end
      %
   end
   %
end

if  strcmp(model.kern.type,'rbfardjit')
    model.learnSigmaf = 1;
end

% initialize inducing inputs by kmeans 
%kmeansops = foptions;
%kmeansops(14) = 10;
%kmeansops(5) = 1;
%kmeansops(1) = 0;
%ch = randperm(model.N);
%centres = model.vardist.means(ch(1:model.k),:);
%model.X_u = kmeans(centres, model.vardist.means, kmeansops);

% Initialise beta according to the desired SNR, which is computed as
% variance_of_centered_data * beta. The relative values of the data
% variance and the beta show how well the model fits the data or explains
% it with noise. The closer 1/beta is to the variance the more information
% the model attempts to explain with noise.
if ~isempty(options) && isfield(options, 'initSNR')
    %if isfield(model, 'mOrig')
    %    mm = model.mOrig;
    %else
    %    mm = model.m;
    %end
    mm = Y;
    if var(mm(:)) < 1e-8
        warning(['Variance in data was too small. Setting beta to 1e+7'])
        model.beta = 1e+7;
    else
        model.beta = 1/((1/options.initSNR * var(mm(:))));
    end
else
    model.beta = 1000;%/max(var(Y));
end


initParams = modelExtractParam(model); %vargplvmExtractParam(model);
model.numParams = length(initParams);
% This forces kernel computation.
model = modelExpandParam(model, initParams); %model = vargplvmExpandParam(model, initParams);

