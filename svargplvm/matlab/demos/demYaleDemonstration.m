% DEMYALEDEMONSTRATION Demonstrating a trained svargplvm model on the Yale faces data, graphically.
%
% SEEALSO: MRDtutorial.m
% COPYRIGHT: Andreas C. Damianou, 2013
% 
% SHEFFIELDML

disp('#--- YALE FACES DEMO ---------')
disp('# Change modality (set of 3 faces) from the variable ''faceSet'' of the demo.')
disp('# For 1st modality, check sampling from dimensions (1,3) and (5,14)')
disp('# For 2nd modality, check sampling from dimensions (1,3) and (10,13)')

%
% The only parameter of the demo: it can take values 1 or 2, depending on
% if we want to visualise the 1st set of 3 faces or the second set of
% another 3 faces (try both!). Each set is a different modality and you can see
% the learned lengthscales in the figure that pops up.

% * Check sampling from shared dimensions 1 and 3 (first set dim. 2 to be
% somewhere on a red cross)
% * For faceSet 1 check sampling from priv. dimensions 5, 14 (largest
% scales). You can fix to a face and return to dims 1,3 to change light!
% * For faceSet 2 check sampling e.g.10 and 13
faceSet = 1;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


experimentNo = 25;

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

if ~exist('itNo')         ,  itNo = [500 1500 1500 1500];              end     % Default: 2000
if ~exist('indPoints')    ,  indPoints = 120;          end     % Default: 49
if ~exist('initVardistIters'), initVardistIters = 180;      end
if ~exist('mappingKern')   ,  mappingKern = {'rbfard2', 'white'}; end

% 0.1 gives around 0.5 init.covars. 1.3 biases towards 0.
if ~exist('vardistCovarsMult'),  vardistCovarsMult=1.3;                  end
% Set to empty value {} to work with toy data
if ~exist('dataSetNames')    ,    dataSetNames = {};    end
if ~exist('invWidthMult'),       invWidthMult = 5;                     end
if ~exist('dataType'), dataType = 'default'; end
if ~exist('latentDimPerModel'), latentDimPerModel = 10; end
if ~exist('experimentNo'), experimentNo = 404; end
if ~exist('doPredictions'), doPredictions = false; end
% If this is true, then the model is in "D > N" mode.
if ~exist('DgtN'), DgtN = false; end
% Create initial X by doing e.g. ppca in the concatenated model.m's or by
% doing ppca in the model.m's separately and concatenate afterwards?
if ~exist('initial_X'), initial_X = 'separately'; end % Other options: 'together'
% Which indices to use for training, rest for test
if ~exist('indTr'), indTr = -1; end

enableParallelism = 0;


% Shuffle one of the two datasets but maintain the correct correspondance,
% i.e. the corresopnding (y_i,z_i) pairs should still be from the same
% angle.
dataType = 'Yale6Sets';
dataSetNames = 'YaleSubset6_1';
[Y,lbls]=svargplvmLoadData(dataSetNames);


N1 = size(Y{1},1);
Yall{2} = [Y{4};Y{5};Y{6}];
identities{2}=[ones(N1,1) 2*ones(N1,1) 3*ones(N1,1)];

numSubsets = 3;
Yall{1} = zeros(numSubsets*N1, size(Y{1},2));
for i=1:N1
    perm = randperm(numSubsets);
    counter = 0;
    for j=perm
        Yall{1}(i+counter*N1,:) = Y{j}(i,:);
        identities{1}(i+counter*N1) = j;
        counter = counter+1;
    end
end


clear Y;


numberOfDatasets = length(Yall);
height = lbls(1); width = lbls(2);


%{
for i=1:size(Yall{1},1)
    for d=1:numberOfDatasets
        subplot(1,numberOfDatasets,d)
        imagesc(reshape(Yall{d}(i,:),height, width)), title(num2str(identities{d}(i))),colormap('gray');
    end
    pause
end
%}




%-- Load datasets
for i=1:numberOfDatasets
    Y = Yall{i};
    dims{i} = size(Y,2);
    N{i} = size(Y,1);
    if indTr == -1
        indTr = 1:N{i};
    end
    indTs = setdiff(1:size(Y,1), indTr);
    Ytr{i} = Y(indTr,:);
    Yts{i} = Y(indTs,:);
    t{i} = linspace(0, 2*pi, size(Y, 1)+1)'; t{i} = t{i}(1:end-1, 1);
    timeStampsTraining{i} = t{i}(indTr,1); %timeStampsTest = t(indTs,1);
    d{i} = size(Ytr{i}, 2);
end

for i=2:numberOfDatasets
    if N{i} ~= N{i-1}
        error('The number of observations in each dataset must be the same!');
    end
end


%%
% demYaleFace03Vargplvm102.mat % --> for one face only (in vargplvm)
load matFiles/demYale6SetsSvargplvm25
model = svargplvmRestorePrunedModel(prunedModel, Ytr);
model.comp{1}.m = model.comp{1}.mOrig;
model.comp{2}.m = model.comp{2}.mOrig;

% bar(model.comp{1}.kern.inputScales)
%set(gca,'YTick',[]);
%set(gca,'FontSize',17)

%%
v = faceSet; % This can be either 1 or 2
lvmVisualiseGeneral(model.comp{v}, [], 'imageMRDVisualise', 'imageMRDModify', false, [height width], 0,0,1);
%{
    modelVis = model.comp{v};
    %bar(model.comp{v}.kern.comp{1}.inputScales);
    figure
    % The following causes OUTOFMEMORY exception except from when we prune the
    % video dimensions: (Note: also, lvmVisualise was changed a bit so that no
    % dynamic slides is presented, because otherwise a strange error occurs).
    modelVis.y = Ytr{v};
    sc = 1; % Set to 4 and try again, if you run out of memory
    [modelP, newHeight, newWidth] = vargplvmReduceVidModel(modelVis, height, width, sc,sc);
    lvmVisualise(modelP, [], 'imageMRDVisualise', 'imageMRDModify', [newHeight newWidth],0,0,1);
    clear modelP
%}
%figure,bar(model.comp{1}.kern.comp{1}.inputScales), figure, bar(model.comp{2}.kern.comp{1}.inputScales)
figure
svargplvmShowScales(model);

fprintf('\n  (You should now have 3 figures openned...)\n')
fprintf('\n# Press any key to continue to solving the correspondence problem...')
pause

%%

%---------------------------- PREDICTIONS ---------------
fprintf(['\n\n # Now we can solve the correspondence problem. Given a test image y*\n', ...
'belonging in one of the modalities'' space, we can find the corresponding\n',...
'latent point x* which the model already has learned how to segment into\n',...
'x* = [x*_y x*_{shared} x*_z], so that:\n',...
'a) if y* actually comes from the training set, we find the training latent\n',...
'point x*_{NN} which is closest to x* in the dimensions x*_{shared} and \n',...
'take the dimensions of x*_{NN} corresponding to the z modality so that to\n',...
'recover the training z*. \n',...
'b) if y* comes from a test set, x* is found by optimisation and z* is\n',...
'generated (novel image) accordingly\n\n']);
reply = input('# Give the number of points to use (0 means skip this part) [10]:', 's');
if isempty(reply)
    reply = '10';
end
numberTestPoints = str2num(reply);

if numberTestPoints ~= 0
    obsMod = 1; % one of the involved sub-models (the one for which we have the data)
    infMod = setdiff(1:2, obsMod);
    
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
    
    if ~exist('testOnTraining')
        testOnTraining=1;
    end
    
    if testOnTraining
        perm = randperm(model.N);
        testInd = perm(1:numberTestPoints);
    else
        perm = randperm(size(Yts{obsMod},1));
        testInd = perm(1:numberTestPoints);
    end
    
    scrsz = get(0,'ScreenSize');
    
    for i=1:length(testInd)
        curInd = testInd(i);
        fprintf('# Testing indice number %d ', curInd);
        if testOnTraining
            fprintf('taken from the training set\n');
            y_star = model.comp{obsMod}.y(curInd,:);
            x_star = model.comp{obsMod}.vardist.means(curInd,:);
            varx_star = model.comp{obsMod}.vardist.covars(curInd,:);
        else
            fprintf('taken from the test set\n');
            y_star = Yts{obsMod}(curInd,:);
            z_star = Yts{infMod}(curInd,:);
            dst = dist2(y_star, model.comp{obsMod}.y);
            [mind, mini] = min(dst);
            
            Init(i,:) = model.vardist.means(mini,:);
            vardistx = vardistCreate(model.comp{obsMod}.vardist.means(mini,:), model.q, 'gaussian');
            vardistx.covars = model.comp{obsMod}.vardist.covars(mini,:);
            model.comp{obsMod}.vardistx = vardistx;
            display=1;
            iters = 250;
            % Find p(X_* | Y_*) which is approximated by q(X_*)
            [x_star, varx_star, modelUpdated] = vargplvmOptimisePoint(model.comp{obsMod}, vardistx, y_star, display, iters);%%%
        end
        numberOfNN = 9;
        % Now we selected a datapoint X_* by taking into account only the
        % private dimensions for Y. Now, based on the shared dimensions of
        % that, we select the closest (in a NN manner) X from the training data.
        fprintf('# Finding the %d NN of X_* with the training X based only on the shared dims.\n', numberOfNN);
        [ind, distInd] = nn_class(model.X(:,sharedDims), x_star(:,sharedDims), numberOfNN, 'euclidean');
        
        ZpredMu = zeros(length(ind), size(model.comp{infMod}.y,2));
        ZpredSigma = zeros(length(ind), size(model.comp{infMod}.y,2));
        
        
        % Find p(y_*|x_*) for every x_* found from the NN
        fprintf('# Predicting images from the NN of X_* ');
        for k=1:numberOfNN
            fprintf('.');
            x_cur = model.X(ind(k),:);
            %x_cur(sharedDims) = x_star(sharedDims); %%% OPTIONAL!!!
            %[ZpredMu(k,:), ZpredSigma(k,:)] = vargplvmPosteriorMeanVar(model.comp{infMod}, model.X(ind(k),:));
            ZpredMu(k,:) = vargplvmPosteriorMeanVar(model.comp{infMod}, x_cur);
        end
        fprintf('\n\n');
        
        
        %-- Plots
        % Open a big figure (first 2 args control the position, last 2 control
        % the size)
      %  figure('Position',[scrsz(3)/100.86 scrsz(4)/6.666 scrsz(3)/1.0457 scrsz(4)/1.0682],...
      figure('units','normalized','outerposition',[0 0 1 1], ...
            'Name',['Fig: ' num2str(i) ' (Exp: ' num2str(experimentNo) ')'],'NumberTitle','off')
      numRows = 2;
        
        if testOnTraining
            numCols = ceil((numberOfNN+1)/numRows);
            plotCounter = 1;
        else
            % For the real test image!
            numCols = ceil((numberOfNN+2)/numRows);
            plotCounter = 2;
        end
        subplot(numRows, numCols, 1)
        imagesc(reshape(y_star,height,width)), title(['Original y (image #' num2str(curInd) ')']), colormap('gray')
        
        if ~testOnTraining
            subplot(numRows, numCols, 2)
            imagesc(reshape(z_star,height,width)), title(['Corresponding z (image #' num2str(curInd) ')']), colormap('gray')
        end
        
        for k=1:numberOfNN
            subplot(numRows, numCols, k+plotCounter)
            imagesc(reshape(ZpredMu(k,:), height, width)), title(['NN #' num2str(k)]), colormap('gray');
        end
        
    end
end

