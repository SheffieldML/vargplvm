function [gradient, delta, delta_percent] = gradchek2(w, func, grad, varargin)
%GRADCHEK Checks a user-defined gradient function using finite differences.
%
%	Description
%	This function is intended as a utility for other netlab functions
%	(particularly optimisation functions) to use.  It enables the user to
%	check whether a gradient calculation has been correctly implmented
%	for a given function. GRADCHEK(W, FUNC, GRAD) checks how accurate the
%	gradient  GRAD of a function FUNC is at a parameter vector X.   A
%	central difference formula with step size 1.0e-6 is used, and the
%	results for both gradient function and finite difference
%	approximation are printed. The optional return value GRADIENT is the
%	gradient calculated using the function GRAD and the return value
%	DELTA is the difference between the functional and finite difference
%	methods of calculating the graident.
%
%	GRADCHEK(X, FUNC, GRAD, P1, P2, ...) allows additional arguments to
%	be passed to FUNC and GRAD.
%   The first additional argument must be the model struct, the second is 
%   optional and must contain the names of the parameters
%
%	See also
%	CONJGRAD, GRADDESC, HMC, OLGD, QUASINEW, SCG
%

%	Copyright (c) Ian T Nabney (1996-2001)

% varargin inputs should be model and names of parameters
assert(length(varargin)==2);

func = fcnchk(func, length(varargin));
grad = fcnchk(grad, length(varargin));

% Treat
nparams = length(w);
step    = zeros(1,length(w));
deltaf = zeros(1, nparams);
epsilon = get_fd_epsilon(w, varargin{2});
for i = 1:nparams
  % Move a small way in the ith coordinate of w
  if( mod(i,100) == 0 )
     fprintf('Processing dimension %d/%d...\n', i, nparams); 
  end
  step(i) = 1.0;
  fplus  = feval('linef', epsilon(i), func, w, step, varargin{1});
  fminus = feval('linef', -epsilon(i), func, w, step, varargin{1});
  % Use central difference formula for approximation
  deltaf(i) = 0.5*(fplus - fminus)/epsilon(i);
  step(i) = 0.0;
end
gradient = feval(grad, w, varargin{1});
fprintf(1, 'Checking gradient ...\n\n');
delta = gradient - deltaf;
delta_percent = (abs(delta) ./ abs(gradient)) * 100;

if numel(varargin) == 2
    names = varargin{2};
    fprintf('\n%5s  %30s  %13s  %13s  %13s  %13s   %2s\n\n', 'n', 'name', 'pval', 'analytic', 'diffs', 'delta', 'percent');
    for i = 1:nparams
        fprintf('%5d  %30s  %13.3e  %13.3e  %13.3e  %13.3e  %2.1f\n', i, names{i}, w(i), gradient(i), deltaf(i), delta(i), delta_percent(i));
    end
else
    fprintf('\n%5s  %13s  %13s  %13s  %13s   %2s\n\n', 'n', 'pval', 'analytic', 'diffs', 'delta', 'percent');
    for i = 1:nparams
        fprintf('%5d  %13.3e  %13.3e  %13.3e  %13.3e  %2.1f\n', i, w(i), gradient(i), deltaf(i), delta(i), delta_percent(i));
    end  
end
fprintf('\n');


end


function epsilon = get_fd_epsilon(params, names)

    nparams = numel(params);
    epsilon = zeros(1,nparams);
    
    count = 1;
    while ~isempty(names)
        strg = regexp(names{1}, '(^[^\(]+)', 'tokens');
        tkns = cellfun(@(x)regexp(x, '(^[^\(]+)', 'tokens'), names);
        idx  = cellfun(@(x)strcmp(x, strg{:}), tkns);
        names(idx) = [];

        nidx = numel(idx);        
        if nidx ~= numel(epsilon)
            %pad with zeros
            idx = logical([zeros(1,nparams-nidx) idx]);
        end
        epsilon(idx) = 1e-3*mean(abs(params(idx)));
        
        count=count+1;
        if count == 50
            error('gradchek2::get_param_groups: Something is wrong with the regex search of the parameter names.');
        end
    end

end

