% UTIL_VARGPLVMCROPVIDEO   
% Note: plot one frame and write "imcrop" on the console. Select the region
% you want to keep and right click -> copy position. this is the cropVector.
% The newHeight and width may be +-1 of
% that, see size(fr) to get it right.

switch dataSetName
        case 'missa'
            % Keep only the head
            newWidth=144;
            newHeight=122;
            Ynew = zeros(size(Y,1), newWidth*newHeight);
            for i=1:size(Y,1)
                fr=reshape(Y(i,:),height,width);
                fr=imcrop(fr,[129.5 91.5 newWidth-1 newHeight-1]);
                Ynew(i,:) = reshape(fr,1,newWidth*newHeight);
            end
            Y = Ynew; width=newWidth; height=newHeight;
        case 'ocean'
            % Keep only the water
            %[1.5 288.5 1279 432]
            newWidth=width;
            newHeight=round(height/1.6666);
            Ynew = zeros(size(Y,1), newWidth*newHeight);
            for i=1:size(Y,1)
                fr=reshape(Y(i,:),height,width);
                fr=imcrop(fr,[1 288.5 newWidth-1 newHeight-1]);
                Ynew(i,:) = reshape(fr,1,newWidth*newHeight);
            end
            Y = Ynew; width=newWidth; height=newHeight;
        case 'horse'
            % Step 1
            cropVector=[24.5 1.5 224 height]; % horse
            %%%
            newHeight=cropVector(4)-1;
            newWidth=cropVector(3)+1;
            Ynew = zeros(size(Y,1), newWidth*newHeight);
            for i=1:size(Y,1)
                fr=reshape(Y(i,:),height,width);
                fr=imcrop(fr,cropVector);
                Ynew(i,:) = reshape(fr,1,newWidth*newHeight);
            end
            Y = Ynew; width=newWidth; height=newHeight;
            
            % Step 2
            cropVector = [0 0 188 159];
            newHeight=cropVector(4);
            newWidth=cropVector(3);
            Ynew = zeros(size(Y,1), newWidth*newHeight);
            for i=1:size(Y,1)
                fr=reshape(Y(i,:),height,width);
                fr=imcrop(fr,cropVector);
                Ynew(i,:) = reshape(fr,1,newWidth*newHeight);
            end
            Y = Ynew; width=newWidth; height=newHeight;
        case 'dog2'
            cropVector=[11.5 3.5 316 357];
            newHeight=cropVector(4);
            newWidth=cropVector(3)+1;
            Ynew = zeros(size(Y,1), newWidth*newHeight);
            for i=1:size(Y,1)
                fr=reshape(Y(i,:),height,width);
                fr=imcrop(fr,cropVector);
                Ynew(i,:) = reshape(fr,1,newWidth*newHeight);
            end
            Y = Ynew; width=newWidth; height=newHeight;
    end  
