dynUsed = 0;



demHighDimVargplvm3

%% ----- Prediction


Nstar = size(YtsOriginal,1);
%reconstrIters = 18000;


predWithMs = 1;

if predWithMs
    %
    display = 1;
    indexP = [];
    Init = [];
    
    fprintf(1, '# Partial reconstruction of test points...\n');
    % w=360; % width
    % h=288; % height
    
    %% OLD
    %{
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
                if exist('cutP')    cutPoint = cutP; end
                mask = [ones(1,cutPoint) zeros(1,h-cutPoint)];
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
    %% %%%%%%%%%%%%%%
    
    
    %}
    
    %%-- NEW:
    demHighDimPrepareTestData
    %%--
    
    mini =[];
    for i=1:size(Yts,1)
        % initialize the latent points using the nearest neighbour
        % from he training data
        dst = dist2(Yts(i,indexPresent), Ytr(:,indexPresent));
        % mind: smaller distance %mini(i): smaller index
        [mind, mini(i)] = min(dst);
    end
    
    
    vardistx = vardistCreate(model.vardist.means(mini,:), model.q, 'gaussian');
    vardistx.covars = model.vardist.covars(mini,:);%0.2*ones(size(vardistx.covars));
    model.vardistx = vardistx;
    
    Yts(:,indexMissing) = NaN;
    iters=reconstrIters;
    
    %-- NEW
    model.mini = mini;
    clear mini %%%
    %---
    [x, varx, modelUpdated] = vargplvmOptimisePoint(model, vardistx, Yts, display, iters);
    
    
    % Find the absolute error
    Testmeans = x;
    Testcovars = varx;
    Varmu = vargplvmPosteriorMeanVar(modelUpdated, x, varx);
    % Mean error per pixel
    errorFull = sum(sum( abs(Varmu(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
    fprintf(1,'# GPLVM Error (in the missing dims) with missing inputs:%d\n', errorFull);
    
    
    
    errorFullPr = sum(sum( abs(Varmu(:,indexPresent) - YtsOriginal(:,indexPresent)) ))/prod(size(YtsOriginal(:,indexPresent)));
    fprintf(1,'# GPLVM Error (in the present dims) with missing inputs:%d\n', errorFullPr);
    
    
    %{
     % Visualization of the reconstruction
         for i=1:Nstar
             subplot(1,2,1);
             fr=reshape(YtsOriginal(i,:),height,width);
             imagesc(fr);
             colormap('gray');
             subplot(1,2,2);
             fr=reshape(Varmu(i,:),height,width);
             imagesc(fr);
             colormap('gray');
             pause(0.5);
         end
    %}
    prunedModelUpdated = vargplvmPruneModel(modelUpdated,1);
    save([fileToSave(1:end-4) 'Pred.mat'], 'Testmeans', 'Testcovars', 'prunedModelUpdated');
    fprintf(1,'# Saved %s\n',[fileToSave(1:end-4) 'Pred.mat']);
    
    
    %-------- NN  ----------
    fprintf(1,'# NN prediction...\n');
    mini =[];
    sortedInd = zeros(size(Yts,1), size(Ytr,1));
    for i=1:size(Yts,1)
        dst = dist2(Yts(i,indexPresent), Ytr(:,indexPresent));
        % mind: smaller distance %mini(i): smaller index
        %[mind, mini(i)] = min(dst);
        [void ,sortedInd(i,:)] = sort(dst,2);
    end
    clear dst %%%
    
    NNmuPart = zeros(Nstar, size(Ytr,2));
    % Set to 1 to find the two NN's in time space. set to 0 to find in
    % data space. timeNN=1 should be set only for datasets created with
    % the "everyTwo" option for the split.
    timeNN=0;
    k=2; % the k parameter in k-NN
    for i=1:Nstar
        if timeNN
            if i < Nstar
                NNmuPart(i, indexMissing) = 0.5*(Ytr(i,indexMissing) + Ytr(i+1,indexMissing));
            else
                NNmuPart(i, indexMissing) = Ytr(end,indexMissing);
            end
        else
            NNmuPart(i,indexMissing) = Ytr(sortedInd(i,1),indexMissing);
            for n=2:k
                NNmuPart(i,indexMissing) = NNmuPart(i,indexMissing)+Ytr(sortedInd(i,n),indexMissing);
            end
            NNmuPart(i,indexMissing) = NNmuPart(i,indexMissing)*(1/k);
        end
        NNmuPart(i, indexPresent) = Yts(i, indexPresent);
    end
    
    % Mean absolute error per pixel
    %errorNNPart = mean(abs(NNmuPart(:) - YtsOriginal(:)));
    errorNNPart = sum(sum( abs(NNmuPart(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
    fprintf(1,'# NN(2) Error (in the missing dims) with missing inputs:%d\n', errorNNPart);
    
    
    %%%%%%  Try different k for k-NN to see which is best
    if ~timeNN
        bestK=2; bestErr = errorNNPart; NNmuPartBest = NNmuPart;
        clear NNmuPart %%%
        for k=[1 3 4 5]; % the k parameter in k-NN, try different values
            NNmuPartK = zeros(Nstar, size(Ytr,2));
            for i=1:Nstar
                NNmuPartK(i,indexMissing) = Ytr(sortedInd(i,1),indexMissing); % first NN
                for n=2:k % more NN's, if k>1
                    NNmuPartK(i,indexMissing) = NNmuPartK(i,indexMissing)+Ytr(sortedInd(i,n),indexMissing);
                end
                NNmuPartK(i,indexMissing) = NNmuPartK(i,indexMissing)*(1/k); % normalise with the number of NN's
                
                NNmuPartK(i, indexPresent) = Yts(i, indexPresent);
            end
            errorNNPartK = sum(sum( abs(NNmuPartK(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
            if errorNNPartK < bestErr
                bestErr = errorNNPartK;
                bestK=k;
                NNmuPartBest = NNmuPartK;
            end
        end
        clear NNmuPartK
    end
    % Mean absolute error per pixel
    fprintf(1,'# NNbest(%d) Error (in the missing dims) with missing inputs:%d\n',bestK, bestErr);
    %%%%%%%%%%%%%
    
    
    %     for i=1:Nstar
    %         subplot(1,2,1);
    %         fr=reshape(YtsOriginal(i,:),height,width);
    %         imagesc(fr);
    %         colormap('gray');
    %         subplot(1,2,2);
    %         fr=reshape(NNmuPart(i,:),height,width);
    %         imagesc(fr);
    %         colormap('gray');
    %         pause(0.5);
    %     end
end