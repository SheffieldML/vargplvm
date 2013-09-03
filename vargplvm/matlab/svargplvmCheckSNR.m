function warnings = svargplvmCheckSNR(SNR, errLimit, warLimit)
% SVARGPLVMCHECKSNR Check Signal to Noise Ratio after
% optimisation, to ensure that the trivial local minimum
% of learning only noise is avoided.
% DESC Check SNR of optimised model
% FORMAT
% ARG SNR: the SNR of optiomised model in a cell array (one cell per
% modality)
% ARG errLimit: Optional, the limit below which an error message
% is printed
% ARG warLimit: Optional, the limit below which a warning message
% is printed
% RETURN warnings: in an array, representing all modalities where
% the SNR is low enough to be considered a warning
%
% COPYRIGHT: Andreas C. Damianou, 2013
%
% VARGPLVM

if nargin < 3 || isempty(warLimit), warLimit = 10; end
if nargin < 2 || isempty(errLimit), errLimit = 2; end
if nargin < 1, error('Not enough arguments given'); end

errStr = sprintf(['\nThis means that a bad local minimum has been reached\n', ...
    'where everything is explained by noise. Please try a different\n', ...
    'initialisation and/or consult the manual.\n']);
warStr = sprintf(['\nThis means that a bad local minimum has been reached\n', ...
    'where everything is explained by noise. Consider trying a different\n', ...
    'initialisation and/or consult the manual.\n']);

errors = [];
warnings = [];
for i = 1:length(SNR)
    if ~isempty(SNR{i}) && SNR{i} <= errLimit
        errors = [errors i];
    end
end

if ~isempty(errors)
    errMsg = 'Error! Low SNR in modalities: ';
    errMsg = [errMsg num2str(errors)];
    errMsg = [errMsg errStr];
    error(errMsg);
else
    for i = 1:length(SNR)
        if ~isempty(SNR{i}) && SNR{i} <= warLimit
            warnings = [warnings i];
        end
    end
end

if ~isempty(warnings)
    warMsg = 'WARNING! Low SNR in modalities: ';
    warMsg = [warMsg num2str(warnings)];
    warMsg = [warMsg warStr];
    warning(warMsg);
end