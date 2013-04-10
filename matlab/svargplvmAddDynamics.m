function  model = svargplvmAddDynamics(model, globalOpt, optionsDyn)
% SVARGPLVMADDDYNAMICS Add a non-standard normal prior to the svargplvm
% model.
% DESC Add a non-standard normal prior to the svargplvm
% model. This prior is referred to as "dynamics", although it can be a
% GP prior of some other type, depending on the inputs to that prior.
% ARG model: The svargplvm model for which the dynamics are added.
% ARG globalOpt: The global options structure with all the demo
% initialisation parameters
% ARG optionsDyn: The options structure for the dynamics field.
% RETURN: the model augmented with a dynamics structure.
%
% SEEALSO : svargplvmCreate
%
% COPYRIGHT: Andreas C. Damianou, 2013
%
% VARGPLVM

    timeConstr = find(strcmp(optionsDyn.constrainType, 'time'));
    dataConstr = find(strcmp(optionsDyn.constrainType, 'data'));
    labelConstr = find(strcmp(optionsDyn.constrainType, 'labels'));
    
    % If the user needs time constraints but didn't provide a time vector,
    % create an equally spaced one.
    if  ~isempty(timeConstr) && (~isfield(optionsDyn, 't') || isempty(optionsDyn.t))
        fprintf(1, '# Time vector unknown; creating random, equally spaced time vector\n');
        %tTemp = linspace(0, 2*pi, size(Y, 1)+1)';
        tTemp = linspace(0, 2*pi, size(model.y, 1)+1)';
        tTemp = tTemp(1:end-1, 1);
        optionsDyn.t = tTemp;
    end
              
    if length(optionsDyn.constrainType) > 1
        if ~strcmp(globalOpt.dynamicKern{1}, 'invcmpnd')
            error('Trying to model more than one constraints with an inappropriate kernel!');
        end        
    end
            
    if ~isempty(dataConstr)
        t{dataConstr} = bc_backConstraintsMatrix(model, globalOpt, options); % model.m; 
    end
    if ~isempty(timeConstr)
        t{timeConstr} = optionsDyn.t;
    end
    
    if ~isempty(labelConstr)
        t{labelConstr} = bc_labelMatrix(optionsDyn.labels, model.y, globalOpt);
    end
    
    % If we have only one constrain then we can pass as an argument to
    % kernCreate a single matrix.
    if length(optionsDyn.constrainType) == 1
        t = t{1};
        optionsDyn.t = t;
    else
        % Although if we have multiple constraints t must be passed as a cell,
        % it then has to be transformed into a concatenated matrix, since the
        % individual constrain kernels will know which part is "theirs".
        optionsDyn.t = [];
        
        for i=1:length(optionsDyn.constrainType)
            optionsDyn.t = [optionsDyn.t t{i}];
        end
    end
    
    %-----------
    
    % Dynamic kernel:
    kern = kernCreate(t, globalOpt.dynamicKern); % Default: {'rbf','white','bias'}

    if globalOpt.initDynKernel
        kern = svargplvmInitDynKernel(kern,globalOpt, optionsDyn);
    end
    
    
    
    optionsDyn.kern = kern;
    
    if isfield(globalOpt,'vardistCovarsMult')
        optionsDyn.vardistCovars = globalOpt.vardistCovarsMult;
    end
    
    % Fill in with default values whatever is not already set
    optionsDyn = vargplvmOptionsDyn(optionsDyn);
    model = vargplvmAddDynamics(model, 'vargpTime', optionsDyn, optionsDyn.t, 0, 0,optionsDyn.seq);
    
    % if not empty, the corresponding gradients of the kernel will be zero,
    % i.e. not learning these elements (typically for the variance of the
    % rbf/matern etc of a invcmpnd)
    if ~isempty(globalOpt.fixedKernVarianceIndices)
        model.dynamics.fixedKernVariance = globalOpt.fixedKernVarianceIndices;
    end
    
    if isfield(globalOpt, 'learnKernelVariance')
        optionsDyn.learnVariance = globalOpt.learnKernelVariance;
    end
    
    
    fprintf(1,'# Further calibration of the initial parameters...\n');
    model = vargplvmInitDynamics(model,optionsDyn);
    
    
    modelAll.optionsDyn = optionsDyn;