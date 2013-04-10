function model = addDefaultVargpTimeDynamics(model, timeStamps, seq)
% ADDDEFAULTVARGPTIMEDYNAMICS
% Depricated. Instead, use vargplvmOptionsDyn. e.g., see
% demStickVargplvmDynMissing1.m
% This file is here for back-compatibility with some older demos.
% VARGPLVM
   
fprintf(1,'# Adding dynamics to the model...\n');

optionsDyn.initX = 'ppca'; %probably not needed

% Create time vector for each dimension; if t is not given, assume an
% equally spaced time vector.
if exist('timeStamps')
    t = timeStamps;
else
    fprintf(1, '# Time vector unknown; creating random, equally spaced time vector\n');
    t = linspace(0, 2*pi, size(model.X, 1)+1)';  
    t = t(1:end-1, 1);
end

%%% Create Kernel -this could go to the vargplvmOptions
kern = kernCreate(t, {'rbf','white'});    

kern.comp{2}.variance = 1e-3; % 1e-1
% The following is related to the expected number of zero-crossings.
kern.comp{1}.inverseWidth = 5./(((max(t)-min(t))).^2);
%kern.comp{1}.inverseWidth = 1;
kern.comp{1}.variance = 1;
optionsDyn.kern = kern;
% The following went to vargpTimeDynamicsCreate.
%optionsDyn.means = model.vardist.means;
%optionsDyn.covars = model.vardist.covars;

if exist('seq')
    model = vargplvmAddDynamics(model, 'vargpTime', optionsDyn, t, 0, 0,seq);
else
    model = vargplvmAddDynamics(model, 'vargpTime', optionsDyn, t, 0, 0);
end
