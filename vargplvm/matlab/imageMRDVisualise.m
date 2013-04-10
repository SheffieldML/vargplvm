function handle = imageMRDVisualise(imageVals, imageSize, transpose, negative, ...
				 scale)

% IMAGEMRDVISUALISE Helper code for showing an image during 2-D visualisation.
% FORMAT
% DESC is a helper function for plotting image data using latent
% variable models. This function is almost identical to imageMorify.m with
% a few additions that allow to save outputs from the GUI for the MRD
% demos.
% ARG imageValues : the values to set the image data to.
% ARG imageSize : the size of the image.
% ARG transpose : whether the resized image needs to be transposed
% (default 1, which is yes).
% ARG negative : whether to display the negative of the image
% (default 0, which is no).
% ARG scale : whether or not to use the imagesc function (defaults
% to 1, which is yes).
% RETURN handle : a the handle to the image data.
%
% COPYRIGHT : Neil D. Lawrence, 2003, 2004, 2006
% MODIFICATIONS: Andreas C. Damianou, 2012
% SEEALSO : imageMRDVisualise, imageModify, lvmResultsDynamic

% VARGPLVM

% Uncomment the following to "binarize" the image (e.g. useful for the
% silhouette demo)

% maxVal = max(max(imageVals));
% minVal = min(min(imageVals));
% thresh = round((maxVal - minVal) /2);
% imageVals(find(imageVals < thresh)) = minVal;
% imageVals(find(imageVals >=thresh)) = maxVal;

if nargin < 3
  transpose = 1;
end
if nargin< 4
  negative = 0;
end
if nargin < 5
  scale = 1;
end
if negative
  imageVals = -imageVals;
end
imageData = reshape(imageVals, imageSize(1), imageSize(2));
if transpose
  imageData = imageData';
end
if scale
  handle = imagesc(imageData);
else
  handle = image(imageData);
end
colormap gray
set(handle, 'erasemode', 'none') % ??

