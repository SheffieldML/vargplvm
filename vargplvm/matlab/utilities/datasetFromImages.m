function [Y, height, width] = datasetFromImages(baseDir, selDirs)

% Create the dataset out of the images.

if nargin < 1
    baseDir=[localDatasetsDirectoryLarge];
end

if nargin < 2
    selDirs = {'normalised_32x32'};
end

for d=1:length(selDirs)
    dirFrom=[baseDir selDirs{d}];
    a=dir(dirFrom);
    counter = 0;
    for i=1:length(a)
        if length(a(i).name)>4 & strcmp(a(i).name(end-2:end),'jpg')
            im = imread([dirFrom filesep a(i).name]);
            %imagesc(im), colormap('gray'); title(a(i).name), pause
            counter = counter+1;
            Yall{d}(counter,:)=im(:)';
        end
    end
    Yall{d} = double(Yall{d});
end
height = size(im,1);
width = size(im,2);

numberOfDatasets = length(Yall);
if d==1
    Y = Yall{1};
end