function handle = imageMRDModify(handle, imageValues, imageSize, transpose, negative, ...
		     scale)

% IMAGEMRDMODIFY Helper code for visualisation of image data.
% FORMAT
% DESC is a helper function for visualising image data using latent
% variable models. This function is almost identical to imageMorify.m with
% a few additions that allow to save outputs from the GUI for the MRD
% demos.
% ARG handle : the handle of the image data.
% ARG imageValues : the values to set the image data to.
% ARG imageSize : the size of the image.
% ARG transpose : whether the resized image needs to be transposed
% (default 1, which is yes).
% ARG negative : whether to display the negative of the image
% (default 0, which is no).
% ARG scale : dummy input, to maintain compatability with
% IMAGEVISUALISE.
% RETURN handle : a the handle to the image data.
%
% COPYRIGHT : Neil D. Lawrence, 2003, 2004, 2006
% MODIFICATIONS: Andreas C. Damianou, 2012
%
% SEEALSO : imageModify, imageMRDVisualise, fgplvmResultsDynamic

% VARGPLVM

%save 'TEMP.mat' 'imageValues'

% Uncomment the following to "binarize" the image (e.g. useful for the
% silhouette demo)
% maxVal = max(max(imageValues));
% minVal = min(min(imageValues));
% thresh = round((maxVal - minVal) /2);
% imageValues(find(imageVals < thresh)) = minVal;
% imageValues(find(imageVals >=thresh)) = maxVal;

if nargin < 4
  transpose = 1;
end
if nargin< 5
  negative = 0;
end
if negative
  imageValues = -imageValues;
end
if transpose
  im = reshape(imageValues, imageSize(1), imageSize(2))';
else
  im = reshape(imageValues, imageSize(1), imageSize(2));
end
set(handle, 'CData', im);

% 
%
if 0
    load TEMPimNumber
    % Initialise TEMPimNumber >100,so that it's alphabetically sorted.
    curNum = num2str(TEMPimNumber);
    tempH = figure;
    imshow(im,[0 180]), colormap('gray') % Def: [0 120];
    axis equal
    axis tight
    axis off %%
    %   text(30,4,'Morphing','HorizontalAlignment','center','VerticalAlignment','top','FontSize',10);
    %'BackgroundColor',[.9 .9 .9], 'FontSize',6);
    %title('Morphing')
    %path = ['../../vargplvm/matlab/Results/diagrams/Yale6Sets/morphingNEWangle/' curNum];
    path = ['/home/andreas/Dropbox/_PhD/Software/github/private/vargplvm/matlab/livMachines/' curNum];
    print( tempH, '-depsc', [path '.eps']) %%%%%%%%%%%
    % print( tempH, '-djpeg100', [path curNum, '.jpeg']) %%%%%%%%%%%
    %imwrite( im, './Results/diagrams/Yale1Face/lightInterp/foo.jpg') %%%%%%%%%%%
                   %  frame = getframe ( gca );
                   % imwrite (frame.cdata,[path '.' 'jpg']); 
    %pause(0.5)
    close(tempH)
    TEMPimNumber = TEMPimNumber +1;
    save('TEMPimNumber.mat', 'TEMPimNumber');
    clear TEMPimNumber
%
end