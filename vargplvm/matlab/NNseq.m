function [NNpredBest, errors] = NNseq(Ytr, YtsOriginal, YtsMissing, windowSizes, startingPoint, replacePresent)

% NNSEQ Sequential NN
% COPYRIGHT: Andreas C. Damianou, 2012
% VARGPLVM

if nargin < 6 || isempty(replacePresent),   replacePresent = true;  end
if nargin < 5 || isempty(startingPoint),    startingPoint = 1;      end
if nargin < 4 || isempty(windowSizes),      windowSizes = 1:5;      end


errors=[];
indexMissing = find(isnan(YtsMissing(1,:)));
indexPresent = setdiff(1:size(YtsMissing,2), indexMissing);



fprintf('  Initialising')
indNN = [];
miniNN = [];
for i=1:size(YtsMissing,1)
    fprintf('.')
    % initialize the latent points using the nearest neighbour from the training data
    dst = dist2(YtsMissing(i,indexPresent), Ytr(:,indexPresent));
    [mind, mini(i)] = min(dst);
    [miniNN(i,:) indNN(i,:)] = sort(dst); %%% mini(i) == indNN(i,1);
end
fprintf('\n')

errInd = 1;
NNpredBest = [];
for ws=windowSizes % Default: 1:5 (ws 11 gives the best results)
    windowSize = ws;
    fprintf('# Testing for windowSize = %d...\n',ws)

    NNpred = Ytr(mini,:);
    NNpredDyn = [];
    mindNN = [];
    miniNN = [];
    dstNN = [];
    NNpredDyn(1,:) = NNpred(startingPoint,:); 
    fprintf('  Predicting')
    % Predict for multiple candidates (as many as the window size) and
    % select the one which is more coherent with the up-to-now sequence.
    for i=2:size(YtsMissing,1)
        fprintf('.')
        candidate = [];
        for j=1:windowSize
            candidate(j,:) = Ytr(indNN(i,j),:);
            % If this is true, then the present dimensions will be replaced
            % by the given ones (as opposed to just returning the whole
            % training NN).
            if replacePresent
                candidate(j,indexPresent) = YtsMissing(i, indexPresent); %% OPTIONAL!! 
            end
        end
        
        if i == 2 %%% OPTIONAL 2 !!!
            % Find the most coherent to the whole first training point (which we
            % know is smooth)
            dstNN = dist2(Ytr(indNN(i,1),:),candidate); %%% OPTIONAL 2 !!!
        else %%% OPTIONAL 2 !!!
            dstNN = dist2(NNpredDyn(i-1,:), candidate);
        end %%% OPTIONAL 2 !!!
        
        [mindNN, miniNN(i)] = min(dstNN);
        NNpredDyn(i,:) = candidate(miniNN(i),:);
    end
   fprintf('\n')
   errors(errInd)= ...
       sum(sum(abs(NNpredDyn(:,indexMissing) - YtsOriginal(:,indexMissing)) ))/prod(size(YtsOriginal(:,indexMissing)));
   
   if errors(errInd) == min(errors)
       NNpredBest = NNpredDyn;
   end
   errInd = errInd+1;
  % meanPose = repmat(mean(Ytr),size(YtsOriginal,1),1);

end
