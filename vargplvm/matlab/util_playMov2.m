function util_playMov2(h, w, options, Y1)


if ~isfield(options,'imFormat')
    imFormat = 'jpg';
end

scrsz = get(0,'ScreenSize');
fntSz=18;

 
        p = options.p;
        if isfield(options, 'indices')
            indices = options.indices;
        end
         if isfield(options, 'sz1') && isfield(options, 'sz2')
             scrsz(3) = options.sz1;
             scrsz(4) = options.sz2;
        end
   
 
    
  
    figure('Position',[scrsz(3)/4.86 scrsz(4)/6.666 scrsz(3)/1.6457 scrsz(4)/1.4682])

 

 if nargin >4 && (length(p) < 2) % Two movies
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
             if isfield(options, 'saveImages')
                 frame = getframe ( gca );
                   imwrite (frame.cdata,[options.saveImages num2str(i) '.' imFormat]); 
                end
        end
    elseif nargin <=4 % One movie
        for i=1:size(Y1,1)
            fr=reshape(Y1(i,:),h,w);
            h=imagesc(fr);
            colormap('gray');
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
            if isempty(p)
                pause
            else
                pause(p);
            end
            if isfield(options, 'saveImages')
               % frame = getframe ( gca );
               %    imwrite (frame.cdata,[options.saveImages num2str(i) '.' imFormat]);
               saveas(h,[options.saveImages num2str(i) '.' imFormat]) 
             end
        end
    else
        fprintf(1,'playMov(h,w,p,Y1,Y2)\n')
    end

