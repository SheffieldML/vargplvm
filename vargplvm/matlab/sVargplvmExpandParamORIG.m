function model = sVargplvmExpandParam(model, params)

%startVal = 1;
%endVal = model.N*model.q;
%paramsX = params(startVal:endVal);
% model.barmu = reshape(params(startVal:endVal), model.N, model.q);

startVal=1;

if isfield(model, 'dynamics') & ~isempty(model.dynamics)
    endVal = model.dynamics.nParams;
else
    endVal = model.vardist.nParams;
end
sharedParams = params(startVal:endVal);



for i = 1:model.numModels
    % model.comp{i}.X = model.X;
    startVal = endVal+1;
    endVal = startVal + model.comp{i}.nPrivateParams-1; 
    params_i = [sharedParams params(startVal:endVal)];
    model.comp{i} = vargplvmExpandParam(model.comp{i}, params_i);
    if strcmp(model.comp{i}.kern.type, 'rbfardjit')
        model.comp{i}.kern.comp{1}.inputScales = model.comp{i}.kern.inputScales;
    end
    %[par, names] = vargplvmExtractParam(model.comp{i}); %%%%5
    %for j=5:-1:1 %%
    %params_i(end-j) %%%
    %end%%
    %model.comp{i}.kern.comp{1}.inputScales
end

% with vargplvmExpandParam for the submodels, and given that the dynamics
% parameters (if there are dynamics) are shared, barmu has been turned into
% mu and this is the same for all submodels. The same holds for the
% variational distribution.
model.X = model.comp{1}.X;

if isfield(model, 'dynamics') & ~isempty(model.dynamics)
    model.dynamics = model.comp{1}.dynamics; % = model.comp{i}.dynamics for all submodels m
end

model.vardist = model.comp{1}.vardist;



