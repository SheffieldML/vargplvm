% DEMHIGHDIMPREPARETESTDATA
% VARGPLVM

    % w=360; % width
    % h=288; % height
    w=width; h=height;
    
    if ~exist('cut')
        cut = 'vertically';
        if strcmp(dataSetName,'missa') || strcmp(dataSetName,'ADN') || strcmp(dataSetName,'claire')
            cut = 'horizontally';
        end
    end
    
    if ~exist('indexMissing')
        switch cut
            case 'horizontally'
                if strcmp(dataSetName,'ADN')
                    cutPoint=round(h/1.55);
                elseif strcmp(dataSetName,'claire')
                    cutPoint=round(h/2);
                elseif strcmp(dataSetName, 'grandma')
                    cutPoint=round(h/1.68);
                else %missa
                    cutPoint=round(h/2)+13;
                end
                if exist('cutP'),    cutPoint = cutP; end
                if ~isscalar(cutPoint)
                    mask = zeros(h,1);
                    mask(cutPoint) = 1;
                else
                    mask = [ones(1,cutPoint) zeros(1,h-cutPoint)];
                    if exist('randomRows') && randomRows
                        perm = randperm(length(mask));
                        mask = mask(perm);
                    end
                end
                mask=repmat(mask, 1,w);
                indexMissing = find(mask);
            case 'vertically'
                if strcmp(dataSetName,'missa')
                    indexMissing=1:round(size(Yts,2)/1.8);
                elseif strcmp(dataSetName,'dog')
                    indexMissing = 1:round(size(Yts,2)/1.70);
                elseif strcmp(dataSetName,'ADN')
                    indexMissing = 1:round(size(Yts,2)/1.7);
                elseif strcmp(dataSetName,'ocean')
                    indexMissing = 1:round(size(Yts,2)/1.6);
                elseif strcmp(dataSetName,'head')
                    indexMissing=1:round(size(Yts,2)/2.08);
                else
                    indexMissing=1:round(size(Yts,2)/2);
                end
                
                % NEW
                if exist('cutP'),    cutPoint = cutP; end
                if ~isscalar(cutPoint)
                    mask = zeros(w,1);
                    mask(cutPoint) = 1;
                else
                    mask = [ones(1,cutPoint) zeros(1,w-cutPoint)];
                    if exist('randomCols') && randomCols
                        perm = randperm(length(mask));
                        mask = mask(perm);
                    end
                end
                mask=repmat(mask, h,1);
                indexMissing = find(mask);
        end
    end
    indexPresent = setdiff(1:model.d, indexMissing);
    
    %
    Yts = YtsOriginal;
    % See the movie after cutting some dims
     %{
     Ytemp=Yts;
     Ytemp(:,indexMissing)=NaN;
     for i=1:size(Ytemp,1)
         fr=reshape(Ytemp(i,:),h,w);
         imagesc(fr);
         colormap('gray');
         pause(0.1)
     end
     clear Ytemp
     %}
    % return %%%%%%%%%%