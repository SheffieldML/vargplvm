% DEMYALEPREDICTONTRAINING Find correspondences between the two subdatasets
% created by the demYaleSvargplvm2 demo.
%
% COPYRIGHT: Andreas C. Damianou, Carl Henrik Ek,  2011
% SEEALSO : demYaleSvargplvm2
%
% VARGPLVM


obsMod = 2; % one of the involved sub-models (possible values: 1 or 2).
infMod = setdiff(1:2, obsMod);


pathToSave = ['/home/andreasd/vargplvm/matlab/Results/CVPR/diagrams/Yale6Sets/grouping/givenMod' num2str(obsMod) '/'];



% Find the dimensions that are shared for obsMod and infMod
if ~exist('sharedDims')
    s1 = model.comp{obsMod}.kern.comp{1}.inputScales;
    s2 = model.comp{infMod}.kern.comp{1}.inputScales;
    % Normalise values between 0 and 1
    s1 = s1 / max(s1);
    s2 = s2 / max(s2);
    
    %  thresh = max(model.comp{obsMod}.kern.comp{1}.inputScales) * 0.001;
    thresh = 0.005;
    
    retainedScales{obsMod} = find(s1 > thresh);
    %thresh = max(model.comp{infMod}.kern.comp{1}.inputScales) * 0.001;
    retainedScales{infMod} = find(s2  > thresh);
    sharedDims = intersect(retainedScales{obsMod}, retainedScales{infMod});
end

% Find X_* only for the shared dimensions (Xs*):
if ~exist('privateDims')
    privateDims = setdiff(1:model.comp{obsMod}.q, sharedDims);
end


% Number of test points to use
numberTestPoints = 60;

perm = randperm(model.N);
testInd = perm(1:numberTestPoints); %%%%%ORIGINAL
%testInd = 160;
testInd = 1:size(Ytr{1},1);

scrsz = get(0,'ScreenSize');

x_star = zeros(length(testInd), size(model.X,2));
for i=1:length(testInd)
    curInd = testInd(i);
    
    y_star = model.comp{obsMod}.y(curInd,:);
    x_star(i,:) = model.comp{obsMod}.vardist.means(curInd,:);
    varx_star = model.comp{obsMod}.vardist.covars(curInd,:);
    
    numberOfNN = 6;
    % Now we selected a datapoint X_* by taking into account only the
    % private dimensions for Y. Now, based on the shared dimensions of
    % that, we select the closest (in a NN manner) X from the training data.
    
    % Actually this is a weighted NN.
    %w = s1(sharedDims)+s2(sharedDims);
    %[ind, distInd] = nn_class2(model.X(:,sharedDims), x_star(i,sharedDims), numberOfNN, 'weighted', w);
    [ind, distInd] = nn_class(model.X(:,sharedDims), x_star(i,sharedDims), numberOfNN, 'euclidean');
    
   % fprintf(['# Given: ' num2str(curInd) '. Found: ']);
   
   %%%% ERRORS
 %{
   fprintf(['Given coords: ' num2str(coordsAll(curInd,:)) '. Found (error): ' ]);
    for k=1:numberOfNN
%         fprintf([num2str(ind(k))]);
%         dif = abs(curInd - ind(k));
%         res = mod(dif, size(Ytr{1}, 1)/3);
%         fprintf(['(' num2str(res) ') ']);

         fprintf([num2str(coordsAll(ind(k),:)) ]);
         
        % fprintf('%.1f ', sqrt(dist2(coordsAll(curInd,:),
        % coordsAll(ind(k),:))));
        
        given = coordsAll(curInd,:);
        found = coordsAll(ind(k),:);
        err(i,k) = abs(given(1) - found(1)) + abs(given(2) - found(2));
        fprintf(['(' num2str(err(i,k)) ') , ']);     
    end
    fprintf('\n');
%}
    %-- Plots
    if exist('makePlots') & ~makePlots
        continue
    end
    
    % Open a big figure (first 2 args control the position, last 2 control
    % the size)
    handle = figure('Position',[scrsz(3)/100.86 scrsz(4)/6.666 scrsz(3)/1.0457 scrsz(4)/1.0682],...
        'Name',['Fig: ' num2str(i) ' (Exp: ' num2str(experimentNo) ')'],'NumberTitle','off');
    numRows = 2;
   
    numCols = ceil((numberOfNN+1)/numRows);
    
    numRows = 1; numCols = 7;%%%%
    plotCounter = 1;
    
    ha = tight_subplot(numRows,numCols,[.01 .01],[.01 .01],[.01 .01]);
    
   % subplot(numRows, numCols, 1)
   axes(ha(1)); 
    im = reshape(y_star, height,width); imshow(im, [0 180]);  
    colormap('gray');% title('Given')
    %imagesc(reshape(y_star,height,width)), title(['Original y (image #' num2str(curInd) ')']), colormap('gray')
    
    
    for k=1:numberOfNN
       % subplot(numRows, numCols, k+plotCounter)
          axes(ha(k+plotCounter)); 
      %  im = reshape(Ytr{infMod}(ind(k),:),height, width); imshow(im, [0  180]); colormap('gray'); % title(['NN #' num2str(k)]) ; %%%
       im = reshape(Ytr{infMod}(ind(k),:),height, width);       
       imshow(im, [0  180]); colormap('gray');  %title(num2str(err(i,k))) ; %%%
        %imagesc(reshape( Ytr{infMod}(ind(k),:), height, width)), title(['NN #' num2str(k)]), colormap('gray');
        axis equal
        axis tight
        axis off
    end
   % print( handle, '-depsc', [pathToSave num2str(curInd) '.eps']) %%%%%%%%%%%
   % system(['epstopdf ' pathToSave num2str(curInd) '.eps']);
    close(handle)
end


