function [model, Ymask]  = plotFmri(prunedModel, dataset, Yplot, options)

% Constants
if nargin < 4
    options = [];
end

if ~isfield(options, 'fixedZdim')
    options.fixedZdim = 15;
end
if ~isfield(options, 'imNo')
    options.imNo = 15;
end




load(dataset);

fName = info.dataset;
dimX = width;
dimY = height;


if strcmp(prunedModel.dataSetInfo.dataSetSplit,'custom')
    Ytr = Y(prunedModel.dataSetInfo.indTr,:);
end
model = vargplvmRestorePrunedModel(prunedModel,Ytr);

% If one output argument exists, return now and don't produce the plots.
if nargout == 1
    return
end

delim = filesep;
p = [localDatasetsDirectoryLarge 'fmri' delim 'fmriDataFinal' delim fName delim 'analyze' delim 'functional' delim 'functional4D.nii'];

if info.applyMask
    maskPath = [localDatasetsDirectoryLarge 'fmri' delim 'fmriDataFinal' delim fName delim 'mask' delim 'lc1ms_deskulled.img'];
    maskNii=load_nii(maskPath);
    Ymask = maskNii.img(:)';
    % If two output arguments exist, return now the model and mask and
    % don't produce any plots.
    if nargout == 2
        return
    end
end


%%
clear maskNii

N = size(Y,1);
nii = load_nii(p);

if nargin < 3 || isempty(Yplot)
    Yplot = Ytr;
end

Yall = zeros(N, size(Ymask,2));
Yall(:,find(Ymask(:)')) = Yplot;
%Yall(:,find(~Ymask(:)')) = 0;

clear Yplot


% Create an alternative dataset with Z direction fixed (so that you can see
% the trajectories over time)
YfixedZ = zeros(N, dimX*dimY);
for i=1:N
    % Store current 3D image by columns
    curImg = nii.img(:,:,:,i);
    % For the alternative dataset
    curImg = nii.img(:,:,options.fixedZdim,i);
    YfixedZ(i,:) = curImg(:)';
end


% Plot image imNo for Z's dimensions dZ=1:end.
% In the other subplot show the variance for ALL images for the
% corresponding Z direction kept fixed
for dZ = 1:dimZ
    subplot(2,1,1)
    resStart = dimX * dimY * (dZ-1) +1;
    resEnd = resStart + dimX * dimY;
    imagesc(reshape(Yall(options.imNo, resStart:resEnd-1), dimX, dimY))
    subplot(2,1,2)
    imagesc(reshape(var(Yall(:,resStart:resEnd-1)),dimX,dimY))
    pause
end

close all
dZ = options.fixedZdim;
for i=1:N
    resStart = dimX * dimY * (dZ)+1;
    resEnd = resStart + dimX * dimY;
    imagesc(reshape(Yall(i, resStart:resEnd-1), dimX, dimY))
    pause
end