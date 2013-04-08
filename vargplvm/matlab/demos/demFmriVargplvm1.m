%... load (from demFmriSvargplvm2.m)
diary LOG_fmri900ProcessedVargplvm222.txt

Y = Yall{1};
keep('Y', 'width', 'height', 'dimZ', 'activity','runNumber');

experimentNo = 222;
dataSetName = 'fmri900Processed'; % 'fmri400Processed3';
mappingKern = 'rbfardjit';
enableParallelism = 0;
indPoints = 120;
latentDim = 60;
dataSetSplit = 'custom';
indTr = 1:size(Y,1);
indTs = size(Y,1);
dynUsed = 0;
initVardistIters = 223;
DgtN = 0;
itNo = [2000 1000 1000 1000 500 500 500 500];
tic;demHighDimVargplvm3;toc
diary off
%%
curMod = 1; % set 1 for the 900 set, 2 for the 400
for i=160:4:250
    x_star = model.X(i,:);
    numberOfNN = 9;
    % Now we selected a datapoint X_* by taking into account only the
    % private dimensions for Y. Now, based on the retained dimensions of
    % that, we select the closest (in a NN manner) X from the training data.
    fprintf('# (%d) Finding the %d NN of X_* with the training X based only on the retained dims.\n', i,numberOfNN);
    [ind, distInd] = nn_class(model.X(:,retainedDims), x_star(retainedDims), numberOfNN, 'euclidean');
    
    fprintf('Real label: %s \nNNs: \n', cell2mat(activity{curMod}(i)));
    for k=1:numberOfNN
        fprintf('  %s\n', cell2mat(activity{curMod}(ind(k))))
    end
    fprintf('\n')
end


