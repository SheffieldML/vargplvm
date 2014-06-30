function handle = imageVisualiseVargplvm(imageVals, imageSize, transpose, negative, ...
				 scale, thresh)

% IMAGEVISUALISE Helper code for showing an image during 2-D visualisation.
% FORMAT
% DESC is a helper function for plotting image data using latent
% variable models.
% ARG imageValues : the values to set the image data to.
% ARG imageSize : the size of the image.
% ARG transpose : whether the resized image needs to be transposed
% (default 1, which is yes).
% ARG negative : whether to display the negative of the image
% (default 0, which is no).
% ARG scale : whether or not to use the imagesc function (defaults
% to 1, which is yes).
% ARG thresh : An array, thresh = [minVal, threshDown, threshUp, maxVal]
% which gives value minVal to all image elements smaller than threshDown
% and value maxVal to all image elements larger than threshUp. 
% RETURN handle : a the handle to the image data.
%
% COPYRIGHT : Neil D. Lawrence, 2003, 2004, 2006
%
% SEEALSO : imageModify, lvmResultsDynamic

% SHEFFIELDML

if nargin < 3 || isempty(transpose)
  transpose = 1;
end
if nargin< 4 || isempty(negative)
  negative = 0;
end
if nargin < 5 || isempty(scale)
  scale = 1;
end
if nargin < 6
    thresh = [];
end
if negative
  imageVals = -imageVals;
end
imageData = reshape(imageVals, imageSize(1), imageSize(2));
if transpose
  imageData = imageData';
end
if ~isempty(thresh)
    imageData(imageData < thresh(2)) = thresh(1);
    imageData(imageData >= thresh(3)) = thresh(4);
end
if scale
  handle = imagesc(imageData);
else
  handle = image(imageData);
end
colormap gray