function channels = demCmu35VargplvmLoadChannels(Ytest, skel)
% DEMCMU35VARGPLVMLOADCHANNELS Given a dataset load the skeleton channels
%
% SEEALSO : demCmu35VargplvmAnimate, demCmu35gplvmVargplvm3.m
% VARGPLVM



% playSkel(Ytest, startInd)
    skel = acclaimReadSkel('35.asf');
    [tmpchan, skel] = acclaimLoadChannels('35_01.amc', skel);
    
    %left indices
    xyzInd = [2];
    xyzDiffInd = [1 3];
    rotInd = [4 6];
    rotDiffInd = [5];
    generalInd = [7:38 41:47 49:50 53:59 61:62];
    startInd = 1;
    endInd = length(generalInd);
    channels(:, generalInd) = 180*Ytest(:, startInd:endInd)/pi;
    startInd = endInd + 1;
    endInd = endInd + length(xyzDiffInd);
    channels(:, xyzDiffInd) = cumsum(Ytest(:, startInd:endInd), 1);
    startInd = endInd + 1;
    endInd = endInd + length(xyzInd);
    channels(:, xyzInd) = Ytest(:, startInd:endInd);
    startInd = endInd + 1;
    endInd = endInd + length(xyzDiffInd);
    channels(:, xyzDiffInd) = cumsum(Ytest(:, startInd:endInd), 1);
    startInd = endInd + 1;
    endInd = endInd + length(rotInd);
    channels(:, rotInd) = asin(Ytest(:, startInd:endInd))*180/pi;
    channels(:, rotInd(end)) = channels(:, rotInd(end))+270;
    startInd = endInd + 1;
    endInd = endInd + length(rotDiffInd);
    channels(:, rotDiffInd) = 0;%cumsum(asin(Ytest(:, startInd:endInd)), 1))*180/pi;
    % skelPlayData(skel, channels, 1/25);
