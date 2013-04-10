% UTIL_SAVEMOVIESCRIPT A demonstration on the use of playMov.m for saving a
% movie.
% DESC Demonstration of how playMov function can be used to save a movie as an
% avi and/or as a set of frames. Also showing how a logo can be placed. The
% script assumes that all the necessary variables are already loaded
% (training set, predictions, height and width of the video).
%
% COPYRIGHT: Andreas C. Damianou, 2011
%
% SEEALSO : util_preprocessVideo.m, util_playMov.m

% VARGPLVM

% The following lines can be commented out if you don't want a logo to be
% added.

 [myLogo,map,alpha]=imread('../logo.png');
options.addLogo = myLogo;
options.logoAlpha = alpha;

% Pause after each frame
options.p = 0.001;

% Assuming that Varmu2 holds the generated images. 
Yall = [Ytr(end-5:end,:);Varmu2];

% Optionally, this will cause the label "Training" or "Generated" to appear
% in each frame
options.indices = zeros(1,size(Yall,1));
options.indices(1:6) = 1;

% To save an 'avi' (better on a windows machine otherwise no compression supported)
% options.saveMovie=1;
options.compression = 'Cinepak';
options.quality = 75;
options.fps = 4;
options.movieName = 'tempMovie.avi';

% To save all frames in a folder
options.saveImages = './temp/dogGenLogo'; % This should be created manually (also: unix('mkdir ./temp'))
options.saveImagesQuality = 100;

playMov(height, width, options, [Ytr(end-5:end,:);Varmu2]);

