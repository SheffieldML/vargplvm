function [Z, newHeight, newWidth] = vargplvmReduceVideo(Yall, h, w, factor1,factor2)

% VARGPLVMREDUCEVIDEO Receives a video and reduces its resolution by factor1 (lines) and factor2).
% DESC Receives a video and reduces its resolution by factor1 (lines) and factor2). The video
% is in a 2-D matrix form.
% FORMAT
% ARG Yall : a matrix for which each line corresponds to a frame. This frame
% must be of dimensions h (height) x w (width) and it must be formed by 
% column-wise serialisation of the original frame, i.e. fr(:)'.
% ARG h : height
% ARG w : height
%
% ARG factor1 : every factor1 lines will be skipped
% ARG factor2: every factor2 columns will be skipped
% RETURN Z : the new video in the same format as described for Yall
% RETURN newHeight: the new height for each frame
% RETURN newWidth: the new width for each frame
%
% COPYRIGHT : Andreas Damianou, 2011
% SEEALSO: vargplvmReduceVidModel.m
% VARGPLVM
  
if nargin == 3
    factor1 = 2;
    factor2 = 2;
elseif nargin == 4
    factor2 = factor1;
elseif nargin > 5 || nargin < 3
    fprintf(1,'# Usage: vargplvmReduceVideo(Yall,h,w,<factor1>,<factor2>\n');
    return
end

 if (factor1 > h-1 || factor2>w-1)
     fprintf(1,'# No changes because factor1 and/or factor 2 > height/width\n');
     return 
 end

for i=1:size(Yall,1)
    Y = Yall(i,:);
    % lines
    zz = zeros(1,factor1);
    mask = [1 zz];
    mask = repmat(mask, 1, ceil(h/(factor1+1)));
    mask = mask(1:h);
    newHeight = sum(mask);
    mask = repmat(mask, 1,w);
    Y = Y(:)';
    Y(find(~mask)) = NaN;

    
    % columns
    zz = zeros(1,h*factor2);
    mask = [ones(1,h) zz];
    mask = repmat(mask, 1, ceil(w/(factor2+1)));
    mask = mask(1:w*h);
    newWidth = floor(sum(mask)/h);
    Y(find(~mask)) = NaN;
    Y = Y(find(~isnan(Y)));
    Z(i,:) = Y;
end

