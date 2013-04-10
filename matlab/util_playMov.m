function util_playMov(h, w, options, Y1, Y2,sz)
% UTIL_PLAYMOV Play a video file which is stored by rows in a 2-D matrix.
% DESC Utility to play video files (a sequence of frames) where each row corresponds
% to one frame, serialized by columns. Each frame must then be reshaped before plotted
% and for this the function must receice its actual dimensions as an argument.
% Sanity warning: This function has been augmented in many steps and is quite messy!!!
% ARG h: the true height of the frame
% ARG w: the true width of the frame
% ARG options: further options for the function
% ARG: Y1 the matrix which keeps the video by rows
% ARG: As Y2, in case the user wants to see two videos side by side
% (as long as both have the same dimensions)
% ARG: sz: The presented frame's dimensions (cf. the true dimensions) is
% divided by this number if the argument is present.
%
% COPYRIGHT: Andreas C. Damianou, 2011
%
% SEEALSO : preprocessVideo.m, saveMovieScript.m

% VARGPLVM


if ~nargin
    fprintf(1,'# Usage: \nplayMov(h,w,[p <p_diff>],Y1,<Y2>,<figSizeFactor>)\n');
    return
end

% The following is only used if options.saveImages exists and is 1.
if ~isfield(options,'imFormat')
    imFormat = 'jpg';
else
    imFormat = options.imFormat;
end

scrsz = get(0,'ScreenSize');
fntSz=18;
 try
    if ~isstruct(options)
        p = options; % pause after each frame
    else
        if ~isfield(options,'p')
            options.p = []; % pause after each frame ([] means that next frame is displayed after pressing any key)
        end
        p = options.p;
        % indices: a vector v in {0,1} with size equal to the number of frames,
        % for which v(i)==1 means it corresponds to a training
        % frame, otherwise to a test or generated frame.
        if isfield(options, 'indices')
            indices = options.indices;
        end
        
        % Save movie as an avi file, also given some more options
        if isfield(options, 'saveMovie') && options.saveMovie
            saveMovie=1;
            if ~isfield(options,'movieName') || ~isfield(options,'fps')
                error('You need to provide the movieName and the fps in the options structure\n');
            end
            aviobj = avifile(options.movieName, 'fps', options.fps); % e.g: 4
            if isfield(options,'quality')
                aviobj.Quality = options.quality;
            end
            if isfield(options,'compression') && ~isempty(options.compression) % e.g. 'Cinepak'
                aviobj.Compression = options.compression;
            end
        end
        
        % sz will cause each frame's size to be reduced by that number.
        if isfield(options, 'sz')
             scrsz(3) = scrsz(3)/options.sz;
             scrsz(4) = scrsz(4)/options.sz;
             fntSz = fntSz/options.sz;
        end
        
        % Like above (sz) but you can re-size the frame with different
        % scales for the horizontal and vertical axis.
         if isfield(options, 'sz1') && isfield(options, 'sz2')
             scrsz(3) = options.sz1;
             scrsz(4) = options.sz2;
        end
    end
    if exist('sz') && sz
        scrsz(3) = scrsz(3)/sz;
        scrsz(4) = scrsz(4)/sz;
    end
    
    
    figure('Position',[scrsz(3)/4.86 scrsz(4)/6.666 1.2*scrsz(3)/1.6457 0.6*scrsz(4)/1.4682])

    % if sz is zero pause to allow user to fix manually the size
    if exist('sz') && ~sz
        pause
    end

    % Oversampling
    if nargin > 4 && (length(p)==2)
        % if p is 2-D then the first element says how much to pause for the
        % second matrix and the second elem how many frames of the second to play in
        % between from the first matrix
        k=1;
        try
            for i=1:size(Y2,1)
                subplot(1,2,1);
                fr=reshape(Y2(i,:),h,w);
                imagesc(fr);
                colormap('gray');
                title(['frame' num2str(i)])
                text(scrsz(4)/6.666,0,['frame' num2str(i)],'HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor',[.7 .9 .7], 'FontSize',fntSz);
                for j=1:p(2)
                    subplot(1,2,2);
                    fr=reshape(Y1(k,:),h,w);
                    k = k+1;
                    imagesc(fr);
                    colormap('gray');
                    drawnow;
                    pause(p(1))
                end
                if isfield(options, 'saveMovie') && options.saveMovie
                    frame = getframe ( gca );
                    aviobj = addframe ( aviobj, frame );
                end
             if isfield(options, 'saveImages')
                 frame = getframe ( gca );
                   imwrite (frame.cdata,[options.saveImages num2str(i) '.' imFormat]); 
                end
            end
        catch e
            e.message
        end
    elseif nargin >4 && (length(p) < 2) % Two movies
        for i=1:size(Y1,1)
            subplot(1,2,1);
            fr=reshape(Y1(i,:),h,w);
            imagesc(fr);
            colormap('gray');
            subplot(1,2,2);
            fr=reshape(Y2(i,:),h,w);
            imagesc(fr);
            colormap('gray');
            title(['frame' num2str(i)])
           % text(scrsz(4)/6.666,0,['frame' num2str(i)],'HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor',[.7 .9 .7], 'FontSize',fntSz);
            if isempty(p)
                pause
            else
                pause(p);
            end
            if isfield(options, 'saveMovie') && options.saveMovie
                frame = getframe ( gca );
                aviobj = addframe ( aviobj, frame );
            end
             if isfield(options, 'saveImages')
                 frame = getframe ( gca );
                   imwrite (frame.cdata,[options.saveImages num2str(i) '.' imFormat]); 
                end
        end
    elseif nargin <=4 % One movie
        for i=1:size(Y1,1)
            fr=reshape(Y1(i,:),h,w);
            if isfield(options, 'addLogo')
              %  possible arguments for the colormap: 64, 256, etc.
               % h1=imshow(grs2rgb(uint8(fr), gray(64)),'Border','tight');
                 h1=imshow(grs2rgb(uint8(fr), colormap('gray')),'Border','tight'); % Default: gray(64)
            else
                 imagesc(fr);
                 colormap('gray');
            end
            drawnow;
            if exist('indices')
                if indices(i)
                    title(['Frame: ' num2str(i) ' (Training)'], 'Color','b','FontSize',fntSz)
                    text(scrsz(4)/6.666,0,['frame' num2str(i) ' (Training)'],'VerticalAlignment','top','BackgroundColor',[.7 .9 .7], 'FontSize',fntSz);
                else
                    title(['Frame: ' num2str(i) ' (Generated)'], 'Color','r','FontSize',fntSz)
                    text(scrsz(4)/6.666,0,['frame' num2str(i) ' (Generated)'],'VerticalAlignment','top','BackgroundColor',[1 0 0], 'FontSize',fntSz);
                end
            else
                title(['frame' num2str(i)])
               % text(scrsz(4)/6.666,0,['frame' num2str(i)],'HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor',[.7 .9 .7], 'FontSize',fntSz);
            end
            if isfield(options,'addLogo')
               % The following only for png and export_fig
               % set(gcf, 'color', 'none'); % Set the figure background to transparent
               % set(gca, 'color', 'none'); % Set the axes background to transparent
                hold on; h2=imshow(options.addLogo,'Border','tight'); axis off; 
                if isfield(options,'logoAlpha')
                    set(h2,'AlphaData',options.logoAlpha); 
                end
                hold off;
            end
                        
            if isempty(p)
                pause
            else
                pause(p);
            end
            if isfield(options, 'saveMovie') && options.saveMovie
                frame = getframe ( gca );
                aviobj = addframe ( aviobj, frame );
            end
            if isfield(options, 'saveImages')
                frame = getframe ( gca );
                imageNo = num2str(i);
                % Pad the image number with leading zeros so that all
                % images have the same length in their filename (useful for
                % sorting frames).
                lengthAll = length(num2str(size(Y1,1)));
                padding = lengthAll - length(imageNo);
                imageNo = [repmat(['0'],1,padding) imageNo];
             %   if isfield(options, 'addLogo')
                    if ~isfield(options, 'saveImagesQuality')
                        options.saveImagesQuality = 75;
                    end
                    if ~isfield(options, 'saveImagesMagn')
                        options.saveImagesMagn =1;
                    end
                    if strcmp(imFormat,'manual')
                        % pause and wait for the user to save it manually
                        fprintf('User saves: %s \n', [options.saveImages imageNo '.' imFormat]);
                        pause
                    else
                         % Use export_fig which crops the borders from each
                        % frame.
                        % if ~strcmp(imFormat, 'eps')
                            export_fig([options.saveImages imageNo], ['-' imFormat], ['-q' num2str(options.saveImagesQuality)], ...
                                ['-m' num2str(options.saveImagesMagn)]);
                        %else
                        %    print('-depsc', [options.saveImages imageNo '.' imFormat]);
                        %    system(['epstopdf ' options.saveImages imageNo '.' imFormat]);
                        %end
                    end
              %  else
               %    imwrite (frame.cdata,[options.saveImages imageNo '.' imFormat]);
                %end
            end
        end
    else
        fprintf(1,'playMov(h,w,p,Y1,Y2)\n')
    end

 catch e
     e.getReport
     if isfield(options, 'saveMovie') && options.saveMovie
         aviobj = close (aviobj);
     end
 end

if isfield(options, 'saveMovie') && options.saveMovie
    aviobj = close ( aviobj );
end
