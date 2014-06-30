function handle = imageModifyVargplvm(handle, imageValues, imageSize, transpose, negative, ...
		     scale,thresh)

% IMAGEMODIFY Helper code for visualisation of image data.
% FORMAT
% DESC is a helper function for visualising image data using latent
% variable models.
% ARG handle : the handle of the image data.
% ARG imageValues : the values to set the image data to.
% ARG imageSize : the size of the image.
% ARG transpose : whether the resized image needs to be transposed
% (default 1, which is yes).
% ARG negative : whether to display the negative of the image
% (default 0, which is no).
% ARG scale : dummy input, to maintain compatability with
% IMAGEVISUALISE.
% ARG thresh : An array, thresh = [minVal, threshDown, threshUp, maxVal]
% which gives value minVal to all image elements smaller than threshDown
% and value maxVal to all image elements larger than threshUp.
% RETURN handle : a the handle to the image data.
%
% COPYRIGHT : Neil D. Lawrence, 2003, 2004, 2006
%
% SEEALSO : imageVisualise, fgplvmResultsDynamic

% SHEFFIELDML


if nargin < 4 || isempty(transpose)
  transpose = 1;
end
if nargin< 5 || isempty(negative)
  negative = 0;
end
if nargin < 6
    thresh = [];
end
if negative
  imageValues = -imageValues;
end
if ~isempty(thresh)
    imageValues(imageValues < thresh(2)) = thresh(1);
    imageValues(imageValues >= thresh(3)) = thresh(4);
end
if transpose
  set(handle, 'CData', reshape(imageValues(1:imageSize(1)*imageSize(2)), imageSize(1), imageSize(2))');
else
  set(handle, 'CData', reshape(imageValues(1:imageSize(1)*imageSize(2)), imageSize(1), imageSize(2)));
end
