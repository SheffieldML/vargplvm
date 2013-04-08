function skelPlayData2(skelStruct, channels, frameLength, lbls)
% SKELPLAYDATA2 Play skel motion capture data for more than one dataset (for comparison).
%
%	Description:
%
%	SKELPLAYDATA2(SKEL, FRAMELENGTH, CHANNELS, lbls) plays channels from a
%	motion capture skeleton and channels. This function just extends the original one by allowing
% 	more datasets to be presented side by side for comparison
%	 Arguments:
%	  SKEL - the skeleton for the motion.
%	  CHANNELS - the channels for each motion arranged in a cell aray.
%	  FRAMELENGTH - the framelength for the motion.
%
%
%	See also
%	BVHPLAYDATA, ACCLAIMPLAYDATA
%	Copyright (c) 2006 Neil D. Lawrence
% 	Modifications: Andreas C. Damianou 2011
% SHEFFIELDML

if nargin < 3
    frameLength = 1/120;
end
clf



scrsz = get(0,'ScreenSize');

%     scrsz(3) = scrsz(3)/sz;
%     scrsz(4) = scrsz(4)/sz;


figure('Position',[scrsz(3)/4.86 scrsz(4)/6.666 scrsz(3)/1.6457 scrsz(4)/1.4682])


figNo = size(channels,2);

for i=1:figNo
    subplot(1,figNo,i);
    gcas{i} = gca;
    handle{i} = skelVisualise(channels{i}(1, :), skelStruct);
    if exist('lbls')
        title(lbls{i})
    end
end



% Get the limits of the motion.
xlim = get(gca, 'xlim');
minY1 = xlim(1);
maxY1 = xlim(2);
ylim = get(gca, 'ylim');
minY3 = ylim(1);
maxY3 = ylim(2);
zlim = get(gca, 'zlim');
minY2 = zlim(1);
maxY2 = zlim(2);
for i = 1:size(channels{1}, 1)
    Y = skel2xyz(skelStruct, channels{1}(i, :));
    minY1 = min([Y(:, 1); minY1]);
    minY2 = min([Y(:, 2); minY2]);
    minY3 = min([Y(:, 3); minY3]);
    maxY1 = max([Y(:, 1); maxY1]);
    maxY2 = max([Y(:, 2); maxY2]);
    maxY3 = max([Y(:, 3); maxY3]);
end
xlim = [minY1 maxY1];
ylim = [minY3 maxY3];
zlim = [minY2 maxY2];

for i=1:figNo
    set(gcas{i}, 'xlim', xlim, ...
        'ylim', ylim, ...
        'zlim', zlim);
end

% Play the motion
for j = 1:size(channels{i}, 1)
    pause(frameLength)
    for i=1:figNo
        subplot(1,figNo,i)
        skelModify(handle{i}, channels{i}(j, :), skelStruct);
    end
end
