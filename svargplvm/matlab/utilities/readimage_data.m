% for the ETH data
function Y = readimage_data(path,type,mask);

current_path = pwd;

cd(path);
d = dir(path);

Y = [];
for(i = 1:1:length(d))
    if(~d(i).isdir)
        if(d(i).name(end-2:end)==type)
            y = imread(d(i).name);
            if(length(size(y))==3)
                y = rgb2gray(y);
            end
            Y = [Y; y(:)'];
        end
    end
end

Y2 = [];
if(mask)
    cd('maps');
    d = dir;
    for(i = 1:1:length(d))
        if(~d(i).isdir)
            if(d(i).name(end-2:end)==type)
                y = imread(d(i).name);
                if(length(size(y))==3)
                    y = rgb2gray(y);
                end
                y = y(:)';
                y = y>0;
                Y2 = [Y2;y];
            end
        end
    end
    Y = Y.*uint8(Y2);
end

cd(current_path);

return;
