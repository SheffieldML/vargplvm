function util_myFprintf(varargin)

fprintf(varargin{:}); % Print to file
tmp = varargin(2:end);
fprintf(1, tmp{:});      % Print to stdout

