%% Demo for semi-supervised learning
% Copyright: Andreas Damianou 2015





%% id of differerent random trials
allTrials = 1:2;
saveDir = './results/semiSupervised_master/';

%%
for trial = allTrials
    % For each trial, try different number of observed outputs
    for NtrObs = [10:10:160]
        keep('saveDir','allTrials','trial','NtrObs','res');
        % Different trials should have different seeds
        curSeed = trial + 1000;
        randn('seed', curSeed);
        rand('seed', curSeed);
        
        dataSetName = 'oil';
        experimentNo = curSeed;
        % This is for when we bootstrap: we try to populate the output
        % labelled space as follows: if Y_i is in the labelled set,
        % (i.e. labels_i exists), then we find this point's embeding, q(x_i) ~
        % N(mu_i, S_i), from which we sample samplesPerObserved points
        % x_samp_i,k, k=1,..., samplesPerObserved and form the pairs
        % {x_samp_i,k, labels_i}, for all k. Here we try two different
        % numbers for samplesPerObserved (a small and a big number).
        samplesPerObserved3 = 6;
        samplesPerObserved4 = 50;
        
        % Number of iterations for initialization of the model (helps
        % avoiding early local minima)
        initVardistIters = 500;
        % Number of optimization iterations
        itNo = [1200 1000];
        
        % Load data
        [Y, lbls] = lvmLoadData(dataSetName);
        labels = transformLabels(lbls)';
        
        % Use Nts test points randomly from the data set. Use the rest
        % NtrObs as training points.
        Nts = 700;
        perm = randperm(size(Y,1));
        ts = perm(1:Nts); t = Nts;
        trObs = perm(t+1:t+NtrObs); t = t+NtrObs;
        trMis = perm(t+1:end);
        
        %%%%% ATTENTION!!!! For compatibility with the rest of the code,
        %%%%% the notation here differs than from the paper:
        %%%%% Y here is Z in paper
        %%%%% labels here is Y in paper
        % SO:
        % YtrObs: Training points in the labelled set (observed and we also not know their corresponding labelsObs
        % YtrMis: Training points in the unlabelled set (observed, but we do not know their labels labelsMis (we have the below but not show them to the model))
        % Yts:    Test points (we'll infer their corresponding labels labelsTs at test time)
        YtrObs = Y(trObs,:); lblsObs = lbls(trObs, :); labelsObs = labels(trObs,:);
        YtrMis = Y(trMis,:); lblsMis = lbls(trMis, :); labelsMis = labels(trMis,:);
        Yts = Y(ts,:); lblsTs = lbls(ts, :);  labelsTs = labels(ts,:);
        
        % Options for bgplvm
        options = vargplvmOptions('dtcvar');
        options.kern = 'rbfardjit';
        options.numActive = min(25, size(YtrObs,1)); % number of ind. points
        options.latentDim = 8; % number of latent dimensions
        options.optimiser = 'scg2';
        options.initSNR = 100; % Initial signal to noise ratio
        
        %----- BGPLVM embed only from the labelled data (no semi-supervised learning here)
        % X is the embedding of only labelled data
        [X, sigma2, W, ~, modelInitVardist] = vargplvmEmbed(YtrObs, options.latentDim, options, initVardistIters, 0 );
        modelInitVardist.globalOpt.experimentNo = experimentNo;
        modelInitVardist.throwSNRError = false;
        model = vargplvmOptimiseModel(modelInitVardist, 1,1, {[], itNo}, 1);
        SNR = vargplvmShowSNR(model,false);
        [x_star_all, varx_star_all, mini] = vargplvmPredictLatent(model, Yts, [], false, 1000, false);
        
        %----- BGPLVM embed from the labelled and unlabelled data (semi-supervised learning)
        % X2 is a riccher embedding because it is learned from more data
        [X2, sigma22, W2, ~, modelInitVardist2] = vargplvmEmbed([YtrObs; YtrMis], options.latentDim, options, initVardistIters, 0 );
        modelInitVardist2.globalOpt.experimentNo = experimentNo;
        modelInitVardist2.throwSNRError = false;
        model2 = vargplvmOptimiseModel(modelInitVardist2, 1,1, {[], itNo}, 1);
        SNR2 = vargplvmShowSNR(model2,false);
        [x_star_all2, varx_star_all2, mini2] = vargplvmPredictLatent(model2, Yts, [], false, 1000, false);
        
        %----  X3 is an even richer embedding, because it is a superset of X2:
        % it starts as X2, but for every point in X2 which has a known
        % label we sample extra points and assign that label. So, X3 is a
        % larger set.
        % First copy X3 = X2 (later we'll augment it)
        X3 = model2.X(1:length(trObs),:);
        % We'll need the variances too in order to sample
        varX3 = model2.vardist.covars(1:length(trObs),:);
        % Copy the labels
        labels3 = labelsObs;
        % Xnew, labelsNew is initialized in size
        Xnew = nan(size(X3,1)*samplesPerObserved3, size(X3,2));
        labelsNew = nan(size(Xnew,1), size(labelsObs,2));
        k=1;
        % Now for each labelled point in the embedding we sample more
        % points
        for n=1:size(X3,1) % At this point X3 == X2
            for kk = 1:samplesPerObserved3
                Xnew(k,:) = X3(n,:) + randn(size(X3(n,:))).*sqrt(varX3(n,:));
                labelsNew(k,:) = labels3(n,:);
                k = k + 1;
            end
        end
        X3 = [X3; Xnew];
        labels3 = [labels3; labelsNew];
        clear 'Xnew' 'labelsNew';
        
        
        
        %---------- As above but obtain many more samples ----------------
        
        %samplesPerObserved4 = min(samplesPerObserved4, floor(10000/NtrObs));
        X4 = model2.X(1:length(trObs),:);
        varX4 = model2.vardist.covars(1:length(trObs),:);
        labels4 = labelsObs;
        Xnew = nan(size(X4,1)*samplesPerObserved4, size(X4,2));
        labelsNew = nan(size(Xnew,1), size(labelsObs,2));
        k=1;
        for n=1:size(X4,1)
            for kk = 1:samplesPerObserved4
                Xnew(k,:) = X4(n,:) + randn(size(X4(n,:))).*sqrt(varX4(n,:));
                labelsNew(k,:) = labels4(n,:);
                k = k + 1;
            end
        end
        X4 = [X4; Xnew];
        labels4 = [labels4; labelsNew];
        clear 'Xnew' 'labelsNew';
        
        
        %% %%%%%%%%%%%%%%%%%%% CLASSIFICATION STEP: %%%%%%%%%%%%%%%%
        % Now we have all 4 embeddings / models. Recap:
        % - Model 1: Embedding only from the labelled data
        % - model 2: Embedding from the labelled and unlabelled data
        % - model 3: As above but bootstrap
        % - model 4: As above but bootstrap more points
        %
        % Given these embeddings, we'll train classifiers from the embedded
        % point to the label space.
        %
        %--------- Train classifiers from X to Z
        %
        %---- Model 1
        % Training of logistic regression classifier (replace with any
        % classifier you want). Do this for each of the embeddings.
        LogRegError = util_fitLogReg(model.X, labelsObs, x_star_all, lblsTs);
        LogRegError2 = util_fitLogReg(model2.X(1:length(trObs),:), labelsObs, x_star_all2, lblsTs);
        LogRegError3 = util_fitLogReg(X3, labels3, x_star_all2, lblsTs);
        LogRegError4 = util_fitLogReg(X4, labels4, x_star_all2, lblsTs);
        % We can also classify not from the embeddings to the labels, but
        % directly from the points. However, when the dimensionality of
        % these points is too big, the following method will automatially
        % do a PCA embedding. Notice that the directly classifying from the
        % points might be better than going through an embedding, if tjhe
        % points' dimensionality is small or the data is relatively
        % simple-structured.
        LogRegErrorOutputs = util_fitLogReg(YtrObs, labelsObs, Yts, lblsTs);
        
        
        %% -------------- Save results in a structure 'res'
        tmp=['tr' num2str(trial) '_NtrObs' num2str(NtrObs)];
        res.(tmp).seed = curSeed;
        res.(tmp).LogRegError = LogRegError;
        res.(tmp).LogRegError2 = LogRegError2;
        res.(tmp).LogRegError3 = LogRegError3;
        res.(tmp).LogRegError4 = LogRegError4;
        res.(tmp).LogRegErrorOutputs = LogRegErrorOutputs;
        res.(tmp).samplesPerObserved3 = samplesPerObserved3;
        res.(tmp).samplesPerObserved4 = samplesPerObserved4;
        res.(tmp).trObs = trObs;
        res.(tmp).ts = ts;
        res.(tmp).trMis = trMis;
        res.(tmp).W = W;
        res.(tmp).W2 = W2;
        res.(tmp).SNR = SNR;
        res.(tmp).SNR2 = SNR2;
        diary(['LOG_demOilSemiSup4.txt']);
        res.(tmp)
        diary off
        save([saveDir 'demOilSemiSup4.mat'], 'res');
    end
end

%% ------------ PLOTS ------------------
close all
cd(saveDir);
ll = [];
lw = 1.7;
ms = 12;
fs = 40;
errorBars = true;
stdMult = 0.8; % Smaller error bar to be in 0

resAll = [];
doneTrials = [];
try
    load(['demOilSemiSup4.mat']);
    ff = fieldnames(res);
    for jj=1:length(ff)
        if isfield(resAll, ff{jj})
            warning('field already exists')
        end
        resAll.(ff{jj}) = res.(ff{jj});
        try
            doneTrials(end+1) = str2num(ff{jj}(3:find(ff{1}=='_')-1));
        catch ee
            error('ll')
        end
    end
catch e
    warning(['Did not find demOilSemiSup4.mat'])
end
doneTrials = unique(doneTrials);


%%% ATTENTION: Reject trials with low SNR (or repeat them with more
%%% initiVardistIters).
% If embedding results in low SNR (i.e. close to 1), then it has to be rejected. There's no
% automatic check for this now, so you can use the next array manually.
excludeTrials = [];
allTrials = setdiff(doneTrials, excludeTrials);
allNtrObs = [10:10:160];
logRegErrorMatrix = NaN(length(allTrials), length(allNtrObs));
logRegErrorSemiSupMatrix = logRegErrorMatrix;
logRegErrorSemiSupMatrix3 = logRegErrorMatrix;
logRegErrorSemiSupMatrix4 = logRegErrorMatrix;
weightKeptMatrix = logRegErrorMatrix;
weigthKeptSemiSupMatrix = logRegErrorMatrix;
SNRMatrix = logRegErrorMatrix;
SNRSemiSupMatrix = logRegErrorMatrix;

thresh = 0.0002;
%

for i=1:length(allTrials)
    for j=1:length(allNtrObs)
        trial = allTrials(i);
        NtrObs = allNtrObs(j);
        try
            logRegErrorMatrix(i, j) = resAll.(['tr' num2str(trial) '_NtrObs' num2str(NtrObs)]).LogRegError;
            logRegErrorSemiSupMatrix(i, j) = resAll.(['tr' num2str(trial) '_NtrObs' num2str(NtrObs)]).LogRegError2;
            logRegErrorSemiSupMatrix3(i, j) = resAll.(['tr' num2str(trial) '_NtrObs' num2str(NtrObs)]).LogRegError3;
            logRegErrorSemiSupMatrix4(i, j) = resAll.(['tr' num2str(trial) '_NtrObs' num2str(NtrObs)]).LogRegError4;
            weightKeptMatrix(i, j) = sum(resAll.(['tr' num2str(trial) '_NtrObs' num2str(NtrObs)]).W > thresh);
            weigthKeptSemiSupMatrix(i, j) = sum(resAll.(['tr' num2str(trial) '_NtrObs' num2str(NtrObs)]).W2 > thresh);
            SNRMatrix(i, j) = resAll.(['tr' num2str(trial) '_NtrObs' num2str(NtrObs)]).SNR;
            SNRSemiSupMatrix(i, j) = resAll.(['tr' num2str(trial) '_NtrObs' num2str(NtrObs)]).SNR2;
        catch e
            warning(['Could not find tr' num2str(trial) '_NtrObs' num2str(NtrObs)])
        end
    end
end

inds = 1:size(logRegErrorMatrix,2);

meanlogRegErrorMatrix = nanmean(logRegErrorMatrix(:,inds));
meanlogRegErrorSemiSupMatrix = nanmean(logRegErrorSemiSupMatrix(:,inds));
meanlogRegErrorSemiSupMatrix3 = nanmean(logRegErrorSemiSupMatrix3(:,inds));
meanlogRegErrorSemiSupMatrix4 = nanmean(logRegErrorSemiSupMatrix4(:,inds));
if errorBars
    ll(end+1)=errorbar(allNtrObs(inds), meanlogRegErrorMatrix, stdMult*nanstd(logRegErrorMatrix(:,inds)), 'x-', 'LineWidth', lw, 'MarkerSize', ms); hold on;
    ll(end+1)=errorbar(allNtrObs(inds), meanlogRegErrorSemiSupMatrix, stdMult*nanstd(logRegErrorSemiSupMatrix(:,inds)), 'ro-', 'LineWidth', lw, 'MarkerSize', ms);
    ll(end+1)=errorbar(allNtrObs(inds), meanlogRegErrorSemiSupMatrix3, stdMult*nanstd(logRegErrorSemiSupMatrix3(:,inds)), 'g+-', 'LineWidth', lw, 'MarkerSize', ms);
    ll(end+1)=errorbar(allNtrObs(inds), meanlogRegErrorSemiSupMatrix4, stdMult*nanstd(logRegErrorSemiSupMatrix4(:,inds)), 'ks-', 'LineWidth', lw, 'MarkerSize', ms);
    xlabel('# Observed')
    ylabel('# Errors')
    legend('Without Z^u','With Z^u', 'With Samples', 'With More Samples');
    xlim([17,163])
    pp = ylim;
    ylim([-20 pp(2)]);
else
    plot(allNtrObs(inds), meanlogRegErrorMatrix, 'x-'); hold on; title('# Errors')
    plot(allNtrObs(inds), meanlogRegErrorSemiSupMatrix, 'ro-');
    plot(allNtrObs(inds), meanlogRegErrorSemiSupMatrix3, 'g+-');  legend('GP','Semi-supervised GP', 'SS-GP with Samples');
    plot(allNtrObs(inds), meanlogRegErrorSemiSupMatrix4, 'ks-');  legend('GP','Semi-supervised GP', 'SS-GP with More Samples');
end

% figure
% plot(allNtrObs(inds), nanmean(weightKeptMatrix(:,inds)), 'x-'); hold on; title(['Weights kept (threshold=' num2str(thresh) ')'])
% plot(allNtrObs(inds), nanmean(weigthKeptSemiSupMatrix(:,inds)), 'ro-');  legend('GP','Semi-supervised GP');
%
% figure
% plot(allNtrObs(inds), nanmean(SNRMatrix(:,inds)), 'x-'); hold on; title('SNR')
% plot(allNtrObs(inds), nanmean(SNRSemiSupMatrix(:,inds)), 'ro-');  legend('GP','Semi-supervised GP');
%



