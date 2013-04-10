function model = vargplvmModelInit(model, globalOpt)

% VARGPLVMMODELINIT Initialise a vargplvm model given global demo options
% COPYRIGHT: Andreas C. Damianou, 2012
% SEEALSO: vargplvmCreate
% VARGPLVM

if globalOpt.DgtN
    model.mOrig = model.m;
    
    model = vargplvmParamInit(model, model.mOrig, model.X);
else
    model = vargplvmParamInit(model, model.m, model.X);
end

% Lengthscales
if strcmp(model.kern.type, 'rbfardjit')
    model.kern.inputScales = globalOpt.invWidthMult./(((max(model.X)-min(model.X))).^2); % Default 5
else
    model.kern.comp{1}.inputScales = globalOpt.invWidthMult./(((max(model.X)-min(model.X))).^2); % Default 5
end



model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));
%model.kern.comp{1}.variance = max(var(Y)); %%%



if model.N > 50 && globalOpt.enableParallelism
    fprintf('# Parallel computations w.r.t the datapoints!\n');
    model.vardist.parallel = 1;
end

%%%
if ~isfield(globalOpt, 'beta') || isempty(globalOpt.betaInit)
    if model.DgtN
        model.beta = 1/((1/globalOpt.initSNR * var(model.mOrig(:))));
    else
        model.beta = 1/((1/globalOpt.initSNR * var(model.m(:))));
    end
else
    model.beta = globalOpt.betaInit;
end

%%% NEW
if model.beta < 1e-7
    warning('Beta was too small... Setting beta to 1e-7')
    model.beta = 1e-7;
elseif model.beta > 1e+7
    warning('Beta was too big... Setting beta to 1e+7')
    model.beta = 1e+7;
end

model.dataSetInfo.dataSetName = globalOpt.dataSetName;

params = vargplvmExtractParam(model);
model = vargplvmExpandParam(model, params);

model.date = date;
model.globalOpt = globalOpt;
