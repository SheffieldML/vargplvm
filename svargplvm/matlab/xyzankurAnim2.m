function xyzankurAnim2(X, fps)

% XYZANKURANIM2 Animate point cloud of stick man from Agarwal & Triggs dataset.
% FORMAT
% DESC animates a matrix of x,y,z point clound positions representing the
% motion of the figure used to generate the silhouttes for Agarwal &
% Triggs silhouette data.
% ARG y : the data to animate.
% ARG fps : the number of frames per second to animate (defaults to 24).
%
% SEEALSO : xyzankurVisualise, xyzankurModify
%  
% COPYRIGHT : Carl Henrik Ek and Neil Lawrence, 2008
% Modifications: Andreas C. Damianou, 2012
%
% Modified from the MOCAP toolbox.
%
% SHEFFIELDML


if(nargin<2)
  fps = 24;
    if(nargin<1)
      error('Too few arguments');
    end
end

for(i = 1:1:size(X,1))
  if(i==1)
    handle = xyzankurVisualise(X(i,:),1);
  else
    xyzankurModify(handle,X(i,:));
  end
  if fps < 0
      pause
  else
    pause(1/fps);
  end
  title(num2str(i))
  %i
end