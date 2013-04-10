function [Y, m, sd] = utils_normaliseData( Y, mn, sdv )

% UTILS_NORMALISEDATA 
% Normalise each dimension to zero mean. This is necessary since the Gaussian processes have zero mean.
% Enforce unit standard deviation for all dimensions jointly
% VARGPLVM

    if nargin == 1
        m  = mean(Y);
        Y  = Y - repmat(m, size(Y,1), 1);
        sd = std(Y(:));
        Y  = Y ./ sd;
    else
        if nargin ~= 3
            error('Invalid number of input arguments.');
        end
        Y = (Y - repmat(mn, size(Y,1), 1)) ./ sdv;
    end

end