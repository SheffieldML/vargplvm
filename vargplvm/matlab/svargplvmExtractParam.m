function [params, names] = svargplvmExtractParam(model)

% SVARGPLVMEXTRACTPARAM Extract a parameter vector from a variational shared GP-LVM model.
% FORMAT
% DESC extracts a parameter vector from a given SVARGPLVM model.
% ARG model : svargplvm model from which to extract parameters
% RETURN params : model parameter vector
%
%
% Parameters are extracted as a vector in the following order (left to right) 
% - parameter{size} -
% vardistParams{model.vardist.nParams} % mu, S
%       OR
% [dynamicsVardistParams{dynamics.vardist.nParams} dynamics.kernParams{dynamics.kern.nParams}] % mu_bar, lambda
%    % Followed by:
% private params of 1st model {model.comp{1}.nPrivateParams}
% private params of 2nd model {model.comp{2}.nPrivateParams}
%           ...
% private params of i-th model {model.comp{i}.nPrivateParams}
%
%
% SEEALSO : svargplvmCreate, svargplvmExpandParam, modelExtractParam
%
% COPYRIGHT : Andreas C. Damianou, 2011

% VARGPLVM

if nargout > 1
    returnNames = true;
else
    returnNames = false;
end

try
    pool_open = matlabpool('size')>0;
catch e
    pool_open = 0;
end

if pool_open && (isfield(model,'parallel') && model.parallel)
    if returnNames
        [params, names] = sVargplvmExtractParamPar(model,returnNames);
    else
        params = sVargplvmExtractParamPar(model,returnNames);
    end
else
    if returnNames
        [params, names] = sVargplvmExtractParamOrig(model,returnNames);
    else
        params = sVargplvmExtractParamOrig(model,returnNames);
    end
end


function [params, names] = sVargplvmExtractParamPar(model, returnNames)
% params = model.X(:)';
if isfield(model, 'dynamics') & ~isempty(model.dynamics)
    if returnNames
        [params, names] = modelExtractParam(model.dynamics);
    else
        params = modelExtractParam(model.dynamics);
    end
else
    if returnNames
        [params, names] = modelExtractParam(model.vardist);
    else
        params = modelExtractParam(model.vardist);
    end
end


% Now extract the "private" parameters of every sub-model. This is done by
% just calling vargplvmExtractParam and then ignoring the parameter indices
% that are shared for all models (because we want these parameters to be
% included only once).
parfor i = 1:model.numModels
    if returnNames
        [params_i,names_i] = vargplvmExtractParam(model.comp{i});
    else
        params_i = vargplvmExtractParam(model.comp{i});
    end
    % params_i = params_i((model.comp{i}.N*model.comp{i}.q)+1:end);
    if isfield(model, 'dynamics') & ~isempty(model.dynamics)
        params_i = params_i(model.dynamics.nParams+1:end);
    else
        params_i = params_i(model.vardist.nParams+1:end);
    end
    paramsC{i} = params_i;
    %params = [params params_i];
    
    if returnNames
        if isfield(model, 'dynamics') & ~isempty(model.dynamics)
            names_i = names_i(model.dynamics.nParams+1:end);
        else
            names_i = names_i(model.vardist.nParams+1:end);
        end
        % names = [names names_i];
        namesC{i} = names_i;
    end
end


for i=1:model.numModels
    params = [params paramsC{i}];
    if returnNames
        names = [names namesC{i}];
    end
end




function [params, names] = sVargplvmExtractParamOrig(model, returnNames)

% params = model.X(:)';
if isfield(model, 'dynamics') & ~isempty(model.dynamics)
    if returnNames
        [params, names] = modelExtractParam(model.dynamics);
    else
        params = modelExtractParam(model.dynamics);
    end
else
    if returnNames
        [params, names] = modelExtractParam(model.vardist);
    else
        params = modelExtractParam(model.vardist);
    end
end

% if returnNames
%     for i=1:size(model.X,1)
%         for j=1:size(model.X,2)
%             if isfield(model, 'dynamics') & ~isempty(model.dynamics)
%                 names{i,j} = ['mubar(' num2str(i) ', ' num2str(j) ')'];
%             else
%                 names{i,j} = ['mu(' num2str(i) ', ' num2str(j) ')'];
%             end
%         end
%     end
%     names = {names{:}};
% end


% Now extract the "private" parameters of every sub-model. This is done by
% just calling vargplvmExtractParam and then ignoring the parameter indices
% that are shared for all models (because we want these parameters to be
% included only once).
rmVardist = false;
for i = 1:model.numModels
    if ~isfield(model.comp{i}, 'vardist') % For memory efficiency
        model.comp{i}.vardist = model.vardist;
        rmVardist = true;
    end

    if returnNames
        [params_i,names_i] = vargplvmExtractParam(model.comp{i});
    else
        params_i = vargplvmExtractParam(model.comp{i});
    end
    % params_i = params_i((model.comp{i}.N*model.comp{i}.q)+1:end);
    if isfield(model, 'dynamics') & ~isempty(model.dynamics)
        params_i = params_i(model.dynamics.nParams+1:end);
    else
        params_i = params_i(model.vardist.nParams+1:end);
    end
    params = [params params_i];
    
    if returnNames
        if isfield(model, 'dynamics') & ~isempty(model.dynamics)
            names_i = names_i(model.dynamics.nParams+1:end);
        else
            names_i = names_i(model.vardist.nParams+1:end);
        end
        names = [names names_i];
    end
    
    if rmVardist
        model.comp{i} = rmfield(model.comp{i}, 'vardist');
    end
end
