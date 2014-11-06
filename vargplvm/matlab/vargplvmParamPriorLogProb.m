function f = vargplvmParamPriorLogProb(model)
% VARGPLVMPARAMPRIORLOGPROB Log. probability of any priors put on parameters of
% a vargplvm model
% DESC Log. pdf of any priors put on parameters of
% a vargplvm model
% SEEALSO: vargplvmAddParamPrior.m, vargplvmParamPriorLogGradients.m,
% priorParamInit.m
% COPYRIGHT: Andreas C. Damianou, 2013
% VARGPLVM

f = 0;
if ~isfield(model, 'paramPriors')
    return
else
    paramPriors = model.paramPriors;
end

for i=1:length(paramPriors)
    cur_params = paramPriors{i}.paramName;
    % cur_ind = paramPriors{i}.index; % That's only for gradients
    switch cur_params
        case 'beta'
            set_param = model.beta;
        case 'input scale'
            if ~isfield(model.kern, 'comp') || strcmp(model.kern.type, 'rbfardjit')
                set_param = model.kern.inputScales;
                %{
                for tr = 1:length(model.kern.transforms)
                    fhandle = str2func([model.kern.transforms(tr).type 'Transform']);
                    if isfield(model.kern.transforms(tr),'transformsettings' ) && ~isempty(model.kern.transforms(tr).transformsettings')
                        set_param = fhandle(model.kern.inputScales, 'xtoa', model.kern.transforms(tr).transformsettings);
                    else
                        set_param = fhandle(model.kern.inputScales, 'xtoa');
                    end
                end
                %}
            else
                set_param = [];
                for cc = 1:length(model.kern.comp)
                    switch model.kern.comp{cc}.type
                        case {'rbfard','rbfard2','rbfardjit'}
                            set_param = [set_param model.kern.comp{cc}.inputScales];
                    end
                    
                end
            end
        case 'dynamicsWhiteKernelVariance'
            for cc = 1:length(model.dynamics.kern.comp)
                if strcmp(model.dynamics.kern.comp{cc}.type, 'white')
                    set_param = model.dynamics.kern.comp{cc}.variance;
                end
            end
    end
    curF = priorLogProb(paramPriors{i}.prior, set_param);
    if isfield(paramPriors{i}.prior, 'scale')
        curF = curF*paramPriors{i}.prior.scale;
    end
    %if isfield(paramPriors{i}.prior, 'log_scale')
    %    curF = curF+paramPriors{i}.prior.log_scale;
    %end
    f = f + curF;
end