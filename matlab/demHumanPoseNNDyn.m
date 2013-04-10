
subsample = 2;
seqToKeep = [1 3 4 5 6];
testSeq = 8;


% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);


if ~exist('indTr'), indTr = -1; end


demHumanPosePrepareData


dataSetNames={'silhouette', 'pose'};
% X_init = Xp;
%mappingKern = {'linard2', 'white'};
%mappingKern = {'rbfard2', 'white'};
latentDim = 5; % Anything > 2 and < 10
%xyzankurAnim(Z_test, 3);
numberOfDatasets = length(Yall);



%-- Load datasets
for i=1:numberOfDatasets
    Y = Yall{i};
    
    dims{i} = size(Y,2);
    N{i} = size(Y,1);
    if indTr == -1
        indTr = 1:N{i};
    end
    %t{i} = linspace(0, 2*pi, size(Yall{i}, 1)+1)'; t{i} = t{i}(1:end-1, 1);
    
    indTs = setdiff(1:size(Y,1), indTr);
    Ytr{i} = Y(indTr,:);
    Yts{i} = Y(indTs,:);
    
    d{i} = size(Ytr{i}, 2);
end
% timeStampsTraining = t{1}(indTr,1); %timeStampsTest = t(indTs,1);


t = linspace(0, 2*pi, length(indTr)+1)'; t = t(1:end-1, 1);

% Fix times:
prevSeq = 1;
timeStampsTraining = [];
dt=0.05;
for i=1:length(seq)
    t = ([0:(seq(i)-prevSeq)].*dt)';
    prevSeq = seq(i)+1;
    timeStampsTraining = [timeStampsTraining ;t];
end;

dt=t(2)-t(1);
timeStampsTest = ([0:size(Y_test,1)-1].*dt)';

%%
errors=[];
for ws=1:5 % Default: 1:5 (ws 11 gives the best results)
    windowSize = ws;
        indNN = [];
        miniNN = [];
       for i=1:size(Z_test,1)
        % initialize the latent points using the nearest neighbour from the training data
        dst = dist2(Y_test(i,:), Ytr{1});
        [mind, mini(i)] = min(dst);
        [miniNN(i,:) indNN(i,:)] = sort(dst); %%% mini(i) == indNN(i,1);
       end
       
    NNpred = Ytr{2}(mini,:);
    NNpredDyn = [];
    mindNN = [];
    miniNN = [];
    dstNN = [];
    initialPoseNN = 1; % Change this to change the initial pose (default: 3)
    NNpredDyn(1,:) = NNpred(initialPoseNN,:); 
    for i=2:size(Z_test,1)
        candidate = [];
        for j=1:windowSize
            candidate(j,:) = Ytr{2}(indNN(i,j),:);
        end
        dstNN = dist2(NNpredDyn(i-1,:), candidate);
        [mindNN, miniNN(i)] = min(dstNN);
        NNpredDyn(i,:) = candidate(miniNN(i),:);
    end
    
    
    meanPose = repmat(mean(Ytr{2}),size(Z_test,1),1);
    errors(ws).meanPose = xyzankurError(meanPose, Z_test);
    errors(ws).NNYspace = xyzankurError(NNpred, Z_test);
    errors(ws).NNDyn = xyzankurError(NNpredDyn, Z_test);
    errorsNNDyn(ws) = errors(ws).NNDyn;
end

%xyzAnkurAnim(NNpredDyn)