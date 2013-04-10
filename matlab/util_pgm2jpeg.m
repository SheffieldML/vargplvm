% PGM2JPEG Converts all pgm files in a directory into jpeg files.
% VARGPLVM

function a = pgm2jpeg(dirFrom, dirTo)
if nargin < 1
    dirFrom = [];
    dirTo = [];
elseif nargin < 2
    dirTo = dirFrom;
end

a=dir(dirFrom);

for i=1:length(a)
    if length(a(i).name)>4 & strcmp(a(i).name(end-2:end),'pgm')
        im = imread([dirFrom filesep a(i).name]);
        curName = [a(i).name(1:end-4) '.jpeg'];
        imwrite(im, [dirFrom filesep curName],'JPEG');
    end
end

%-- Examples
%{
% e.g .1:
a=pgm2jpeg([localDatasetsDirectoryLarge
    'CroppedYale/CroppedYale/yaleB03']);

% e.g. 2:
for j=4:9
    curFolder = ['yaleB0' num2str(j)];
    pgm2jpeg([localDatasetsDirectoryLarge 'CroppedYale/CroppedYale/' curFolder]);
end
for j=10:39
    curFolder = ['yaleB' num2str(j)];
    pgm2jpeg([localDatasetsDirectoryLarge 'CroppedYale/CroppedYale/' curFolder]);
end
%}