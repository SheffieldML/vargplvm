function handle = xyzankurDrawColor(joint, handle, symb, pos, allDirections)

% XYZANKURDRAW Helper function for drawing the point cloud from Agarwal and Triggs data.
% FORMAT
% DESC draws the stick figure associated with a set of x,y,z points for
% the Agarwal and Triggs silhouette data.
% ARG joint : the positions of the joints.
% ARG handle : graphics handles to modify.
% RETURN handle : modified graphics handles.
%
% SEEALSO : xyzankurModify, xyzankurVisualise
%
% COPYRIGHT : Carl Henrik Ek and Neil Lawrence, 2008

% SHEFFIELDML

if nargin < 5 || isempty(allDirections)
    allDirections = false;
end

%zeroinf = 1:6;
%xinf = 7:3:63;
%yinf = 8:3:63;
%zinf = 9:3:63;

limb{1} = [2 3;3 4;4 5]; %spine
limb{2} = [3 6;6 7;7 8;8 9]; % left-arm
limb{3} = [3 10;10 11;11 12;12 13]; % right-arm
limb{4} = [2 14;14 15;15 16;16 20]; % left-leg
limb{5} = [2 17;17 18;18 19;19 21]; % right-leg


if(nargin<2)
  % draw figure
  k = 1;
  for(i = 1:1:length(limb))
    linestyle = '-';
    markersize = 20;
    marker = '.';
    for(j = 1:1:size(limb{i},1))
      if (i==3 && j>3) || (i==5 && j>3)
        markersize = 5;
        marker = 's';
      end
      handle(k) = line(joint(limb{i}(j,:),1),joint(limb{i}(j,:),2), ...
                       joint(limb{i}(j,:),3), 'LineWidth', 3, ...
                       'LineStyle', linestyle, ...
                       'Marker', marker, ...
                       'markersize', markersize);
      k = k + 1;
    end
    end
else
  % modify figure
  k = 1;
  for(i = 1:1:length(limb))
    for(j = 1:1:size(limb{i},1))
        set(handle(k),'Xdata',joint(limb{i}(j,:),1),'Ydata',joint(limb{i}(j,:),2),'Zdata',joint(limb{i}(j,:),3));
        k = k+1;
    end
  end
end
hold on

tt=[0.6 0.6 1.1];
counter = 1;
% modify figure
for k=1:3:61
    if allDirections
        plot3(pos(k+2), -pos(k), -pos(k+1), [symb{counter}], 'MarkerSize',16, 'LineWidth',3);
        counter = counter + 1;
    else
        plot3(pos(k+2)+tt(1), -pos(k)+tt(2), -pos(k+1)+tt(3), [symb{k}(2) 'r'], 'MarkerSize',14, 'LineWidth', 2)
        plot3(pos(k+2)-tt(1), -pos(k)-tt(2), -pos(k+1)-tt(3), [symb{k+1}(2) 'k'], 'MarkerSize',14, 'LineWidth', 2)
        plot3(pos(k+2),       -pos(k),       -pos(k+1),       [symb{k+2}(2) 'm'], 'MarkerSize',14, 'LineWidth', 2)
    end
end

