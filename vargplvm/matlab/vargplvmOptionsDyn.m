function optionsDyn = vargplvmOptionsDyn(optionsDyn, X)
% VARGPLVMOPTIONSDYN Fill in an options structure with default parameters
% FORMAT
% DESC takes an VARGPLVM options structure amends it with default values
% wherever the field is empty.
% ARG optionsDyn : the VARGPLVM structure to fill in with default values.
% ARG X : a latent space initialised for the basic model (only used if
% times t are not given)
%
% COPYRIGHT : Michalis K. Titsias, 2011
% COPYRIGHT : Neil D. Lawrence, 2011
% COPYRIGHT : Andreas C. Damianou, 2011
% 
% SEEALSO : vargplvmInitDynamics, vargplvmAddDynamics

% VARGPLVM


if isfield(optionsDyn, 'type') && ~isempty(optionsDyn.type)
    optionsDyn.type = 'vargpTime';
end


if strcmp(optionsDyn.type, 'vargpTime')
    if ~isfield(optionsDyn, 'initX')
        optionsDyn.initX = 'ppca';
    end
    
    if  ~isfield(optionsDyn, 'constrainType') || isempty(optionsDyn.constrainType)
        % If the dynamics kernel is of type invcmpnd, this cell will
        % explain what type of constraints each element will impose (if
        % other kernel is used, then this cell will only have one element).
        % Possible values (so far) include: "time" and "data".
        optionsDyn.constrainType = {'time'};
    end
    
    % Create time vector for each dimension; if t is not given, assume an
    % equally spaced time vector.
    if  ~isfield(optionsDyn, 't') || isempty(optionsDyn.t)
        fprintf(1, '# Time vector unknown; creating random, equally spaced time vector\n');
        t = linspace(0, 2*pi, size(X, 1)+1)';
        t = t(1:end-1, 1);
        optionsDyn.t = t;
    else
        t=optionsDyn.t;
    end
    
    
    % A better way to initialize the  kernel hyperparameter,
    % especially lengthscale, should be to learn it by running few iterations of
    % GP regression marginal likelihood maximization given as data the PCA output
    % (the default value is jsut reasonable and it sets the inverse lenthscale to quite
    % small value so the GP dynamic prior is weaker (less smoothed) at
    % initializations
    if  ~isfield(optionsDyn, 'kern') || isempty(optionsDyn.kern)
       % kern = kernCreate(t, {'rbf','white'});
        kern = kernCreate(t, {'rbf','white', 'bias'});
        kern.comp{2}.variance = 1e-3; % 1e-1
        
        if  isfield(optionsDyn, 'inverseWidth')
            invWidth=optionsDyn.inverseWidth;
        else
            invWidth=5;
        end
        
        % The following is related to the expected number of zero-crossings.
        kern.comp{1}.inverseWidth = invWidth./(((max(t)-min(t))).^2);
        kern.comp{1}.variance = 1;
        
        optionsDyn.kern = kern;
    end
    
    
    if ~isfield(optionsDyn,'seq')
        optionsDyn.seq=[];
    end
    
    % Set to 1 to reoptimise all inducing points during test time
     if ~isfield(optionsDyn, 'testReoptimise')
     	optionsDyn.testReoptimise = 1;  % !!! SET THE DEFAULT TO 1
     end
    
     % Set to 0 so that if you use a kernel from the exp. family for the
     % dynamics, then its variance is not learned and the lenghtscale
     % compensates this (we have one less free parameter then).
    if ~isfield(optionsDyn, 'learnVariance')
        optionsDyn.learnVariance = 0; % !!! SET THE DEFAULT TO 0
    end
    
    % If set to 1 (not recommended) then the initial means have some
    % smoothing in the first and last values (i.e. closer to zero). Useful
    % mostly for periodic kernels.
    if ~isfield(optionsDyn, 'regularizeMeans')
        optionsDyn.regularizeMeans = 0; 
    end
    
end

