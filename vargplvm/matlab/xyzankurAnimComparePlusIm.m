function xyzankurAnimComparePlusIm(X, X2, Im, sh, fps, names)

% XYZANKURANIMCOMPAREPLUSIM Animate a prediction and ground truth for stick
% man from Agarwal & Triggs dataset and include the corresponding
% silhouette.
% FORMAT
% DESC animates a matrix of x,y,z point clound positions representing the
% motion of the figure used to generate the silhouttes for Agarwal &
% Triggs silhouette data. 
% ARG X : the data to animate.
% ARG X2 : other data to animate.  
% ARG Im : A matrix where each row is a frame of the silhouette modality.
% ARG sh : A vector, such that sh(1) is the height and sh(2) the width of
% the silhouette image (i.e. the number of pixels)
% ARG fps : the number of frames per second to animate (defaults to 24).
% ARG names: names for the modalities in X, X2
%
% SEEALSO : xyzankurVisualise, xyzankurModify
%  
% COPYRIGHT : Carl Henrik Ek and Neil D. Lawrence, 2008
% MODIFICATIONS: Alfredo Kalaitzis, 2012, Andreas Damianou, 2014

% SHEFFIELDML

nDatasets = length(X2);
more_handles = zeros(nDatasets,19);
if isempty(fps), fps=24; end

clf
for i = 1:1:size(X,1)
    if (i == 1)
        subplot(1, nDatasets+2, 1), title(names{1})
        handle = xyzankurVisualise(X(i,:), 1);
        for j = 1:nDatasets
            subplot(1, nDatasets+2, j+1), title(names{j+1})
            more_handles(j,:) = xyzankurVisualise(X2{j}(i,:), 1);
        end
        subplot(1,nDatasets+2, nDatasets+2), title('silhouette')
        handleIm = imagesc(reshape(Im(i,:), sh(1), sh(2))); colormap('gray');
        axis off
        fprintf(1, 'Press ANY key...\n');
        pause
    else
        xyzankurModify(handle, X(i,:));
        for j = 1:nDatasets
            xyzankurModify(more_handles(j,:), X2{j}(i,:));
        end
        set(handleIm, 'CData', (reshape(Im(i,:), sh(1), sh(2))));
    end
    if fps > 0
        pause(1/fps);
    else
        pause
    end
end

% 
% 
%   
%   if(nargin<5)
%     fps = 24;
%     if(nargin<2)
%       fid = 1;
%       if(nargin<1)
%         error('Too few arguments');
%       end
%     end
%   end
%   clf
%   for(i = 1:1:size(X,1))
%     if(i==1)
%       subplot(1, 3, 1)
%       handle = xyzankurVisualise(X(i,:), 1);
%       subplot(1, 3, 2)
%       handle2 = xyzankurVisualise(X2(i,:), 1);
%       subplot(1,3,3)
%       handle3 = imagesc(reshape(Im(i,:), sh(1), sh(2))); colormap('gray');
%     else
%       xyzankurModify(handle, X(i,:));
%       xyzankurModify(handle2, X2(i,:));
%       set(handle3, 'CData', (reshape(Im(i,:), sh(1), sh(2))));
%     end
%     pause(1/fps);
%   end
% end
% 
% 
% 
