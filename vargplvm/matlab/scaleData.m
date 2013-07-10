function [Y, bias, scale] = scaleData(Y, scaleMethod, scaleVal)
% SCALEDATA Scale and center data
% VARGPLVM

if nargin < 3, scaleVal = [];   end
if nargin < 2, scaleMethod = []; end

d = size(Y,2);
bias = mean(Y);

% Remove bias
m = Y;
for i = 1:d
  m(:, i) = m(:, i) - bias(i);
end



scale = ones(1, d);

if ~isempty(scaleMethod) && scaleMethod ~= 0
    if ~isempty(scaleVal) && scaleVal
        warning('Both scale2var1 and scaleVal set for GP');
    end
    if(scaleMethod == 1) % Scale to variance 1
        scale = std(Y);
    elseif scaleMethod == 2 % Scale so that maximum is 1
        scale = max(max(abs(m)));
    else
        error('Unknown scale option')
    end
    scale(find(scale==0)) = 1;
end
if ~isempty(scaleVal) && scaleVal
    scale = repmat(scaleVal, 1, d);
end

% Apply scale
for i = 1:d
  if scale(i)
    m(:, i) = m(:, i)/scale(i);
  end
end

Y = m;