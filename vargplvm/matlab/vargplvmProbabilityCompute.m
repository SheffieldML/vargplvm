function logProb = vargplvmProbabilityCompute(model, y, display, iters)

% VARGPLVMPROBABILITYCOMPUTE description
  
% VARGPLVM
  
% Takes an input a trained vargplvm and a test data point (with possibly missing values)
% Computes the probability density in the test data point 


% Indices of missing dimension
indexMissing = isnan(y);
indexPresent = ~indexMissing;
nTestSamples = size(y,1);
logProb = zeros(nTestSamples,1);

% if ~matlabpool('size')
%     if display, fprintf('vargplvmProbabilityCompute.m: line 17 -- Starting up MatLab pool for parallel computing.'); end
%     matlabpool
% end
try
    nWorkers = matlabpool('size');
catch e
    nWorkers = 1;
end
if nWorkers==0 % If matlabpool has not started, don't try to use it.
    nWorkers = 1;
end

% compute the variational lower without the new data point 
Fold = vargplvmLogLikelihood(model);

if all(all(repmat(indexMissing(1,:), [nTestSamples, 1]) == indexMissing))
    % Faster option:
    if display, 
        fprintf('Optimisation will be done all at once as the missing variables are always the same (or if there is no missing variable).\n');
        fprintf('Initialising the latent point using the nearest neighbour from the training data.\n');
    end
    dst = dist2(y(:,indexPresent(1,:)), model.y(:,indexPresent(1,:)));
    [~, mini] = min(dst,[],2);

    if display, fprintf('Creating the variational distribtion for the test latent point.\n'); end
    vardistx = vardistCreate(model.vardist.means(mini,:), model.q, 'gaussian');

    if display, fprintf('Optimising over the latent points.\n'); end
    model.vardistx = vardistx;
    [X, varX] = vargplvmOptimisePoint(model, vardistx, y, display, iters);

    if display, fprintf('Computing the variational lower with the new data points included.\n'); end
    vardistx.means = X; 
    vardistx.covars = varX;

    if display, fprintf('Generating results for each sample separately (so optimisation is done only once for all test samples).\n'); end
    parfor t = 1:nTestSamples,  
        vardistxTmp = vardistx;
        vardistxTmp.numData = 1;
        % if display, fprintf('Testing sample %d/%d\n', t, nTestSamples);
        vardistxTmp.means = vardistx.means(t,:);
        vardistxTmp.covars = vardistx.covars(t,:) ;
        Fnew = vargplvmPointLogLikelihood(model, vardistxTmp, y(t,:));

        % compute the probability 
        logProb(t) = Fnew - Fold;
    end
    
else % This version is slower, but also uses multithreading:
    c = clock;
    tmpDir = sprintf('/tmp/%s/%d_%d_%d_%d_%d', getenv('USER'), c(1), c(2), c(3), c(4), c(5));
    if display, 
        fprintf('There are missing variables but they are not consistent. Optimising for each test point separately.\n');
        fprintf('In order to reduce the memory usage by each worker, I will use temporary files in %s\n', tmpDir);
    end
    if ~exist(tmpDir,'dir')
        mkdir(tmpDir);
    end

    logProbCells = cell(nWorkers,1);
    for w=1:nWorkers
        tmpFile = sprintf('%s/Ytest_cpu%02d.mat', tmpDir, w);
        if ~exist(tmpFile,'file')
            currSamplesBegin = 1 + (w-1)*floor(nTestSamples/nWorkers);
            currSamplesEnd = w*floor(nTestSamples/nWorkers);
            if w==nWorkers
                currSamplesEnd = nTestSamples;
            end
            sampleIndices = [currSamplesBegin currSamplesEnd];
            if display, fprintf('Writing samples for worker %d/%d: from %d to %d...', w, nWorkers, currSamplesBegin, currSamplesEnd); end
            currY = y(currSamplesBegin:currSamplesEnd,:);
            currIndexPresent = indexPresent(currSamplesBegin:currSamplesEnd,:);
            save(tmpFile, 'currY', 'currIndexPresent', 'sampleIndices');
            if display, fprintf(' done!\n'); end
        end
    end
    clear y;
    clear indexPresent;
    
    parfor w=1:nWorkers,
      tmpFile = sprintf('%s/Ytest_cpu%02d.mat', tmpDir, w);
      if display, fprintf('Loading %s...', tmpFile); end
      testData = load(tmpFile);
      if display, fprintf(' done!\n'); end
      
      currSamples = (testData.sampleIndices(1):testData.sampleIndices(end))';
      currLogProb = zeros(currSamples(end)-currSamples(1)+1,1);

      for tIdx=1:length(currSamples),
        modelCopy = model;
        t = currSamples(tIdx);
        if display, fprintf('-----\nvargplvmProbabilityCompute(): evaluating sample %d of %d\n', t, nTestSamples); end;
        
        dst = dist2(testData.currY(tIdx,testData.currIndexPresent(tIdx,:)), modelCopy.y(:,testData.currIndexPresent(tIdx,:)));
        [~, mini] = min(dst);
        
        % create the variational distribtion for the test latent point
        vardistx = vardistCreate(modelCopy.vardist.means(mini,:), modelCopy.q, 'gaussian');
       
        % optimize over the latent point
        modelCopy.vardistx = vardistx;
        [X, varX] = vargplvmOptimisePoint(modelCopy, vardistx, testData.currY(tIdx,:), display, iters);
            
        % compute the variational lower with the new data point included
        vardistx.means = X; 
        vardistx.covars = varX;

        % Old way:
        Fnew = vargplvmPointLogLikelihood(modelCopy, vardistx, testData.currY(tIdx,:));

        % compute the probability 
        currLogProb(tIdx) = Fnew - Fold;
      end
      logProbCells{w} = currLogProb;
    end
    
    % Concatenating the results of all workers:
    for w=1:nWorkers
      currSamples = ((1 + (w-1)*floor(nTestSamples/nWorkers)) : w*floor(nTestSamples/nWorkers))';
      if w==nWorkers && (currSamples(end) < nTestSamples),
          currSamples = [currSamples; ((currSamples(end)+1):nTestSamples)'];
      end
      logProb(currSamples) = logProbCells{w};
    end
end
end
