function handle = xyzankurVisualise2(pos, v)

% XYZANKURVISUALISE2 Draw the Agarwal & Triggs figure return the graphics handle.
% FORMAT
% DESC draws the stick figure from the Agarwal and Triggs silhouette
% data. 
% ARG pos : the positions of the joints.
% ARG v : the view point for the figure (defaults to standard 3D view). 
% RETURN handle : the graphics handle of the drawn object.
%
% SEEALSO : xyzankurDraw, xyzankurModify
%
% COPYRIGHT : Carl Henrik Ek and Neil Lawrence, 2008
% Modifications: Andreas C. Damianou, 2012
%
% Modified from the MOCAP toolbox.
%
% VARGPLVM



% Convert positions for plotting.
joint = xyzankur2joint(pos);
handle = xyzankurDraw(joint);
axis equal
set(gca,'XLim',[-15 15],'YLim',[-15 15],'ZLim',[0 70]);
if nargin < 2
    view(3)
else
    view(v)
end

if(exist('v','var')&&length(v)==2)
  view(v(1),v(2));
end
xlabel('x')
ylabel('z')
zlabel('y')
axis on

return;