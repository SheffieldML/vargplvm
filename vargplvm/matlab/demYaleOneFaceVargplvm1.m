% DEMYALEONEFACEVARGPLVM1 Run the variational GP-LVM on one of the faces in the YALE dataset.
% COPYRIGHT: Andreas C. Damianou, Carl Henrik Ek, 2011

% VARGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

% Define constants (in a manner that allows other scripts to parametrize
% this one).
if ~exist('experimentNo') ,  experimentNo = 404;      end
if ~exist('itNo')         ,  itNo = [600 600 600];   end     % Default: 2000
if ~exist('indPoints')    ,  indPoints = 64;          end     % Default: 49
if ~exist('latentDim')    ,  latentDim = 7;          end
% Set to 1 to use dynamics or to 0 to use the standard var-GPLVM
if ~exist('dynUsed')      ,  dynUsed = 0;             end
if ~exist('initVardistIters'), initVardistIters = 180;      end     % DEFAULT: 23
if ~exist('dynamicKern')   ,     dynamicKern = {'rbf', 'white', 'bias'}; end
if ~exist('mappingKern')   ,     mappingKern = 'rbfardjit'; end %{'rbfard2', 'bias', 'white'}; end
if ~exist('vardistCovarsMult'),  vardistCovarsMult=1.3;                  end
if ~exist('invWidthMultDyn'),    invWidthMultDyn = 100;                     end
if ~exist('invWidthMult'),       invWidthMult = 5;                     end
if ~exist('initX'),     initX ='ppca';   end
if ~exist('DgtN'), DgtN = 1; end

svargplvmModel = 'demYale6SetsSvargplvm24';
if ~exist('svargplvmYtr')
    %svargplvmYtr = load(''); % ASSUME it's already loaded.
    load svargplvmYtr
end

% Create the dataset out of the images.
baseDir=[localDatasetsDirectoryLarge 'CroppedYale' filesep 'CroppedYale'];
selDir = '03';  %{'03', '08'}; % 03,08

dirFrom=[baseDir filesep 'yaleB' selDir];
a=dir(dirFrom);
counter = 0;
for i=1:length(a)
    if length(a(i).name)>4 & strcmp(a(i).name(end-2:end),'pgm') ...
            & ~strcmp(a(i).name(end-10:end-4),'Ambient')
        im = imread([dirFrom filesep a(i).name]);
        %imagesc(im), colormap('gray'); title(a(i).name), pause
        counter = counter+1;
        Y(counter,:)=im(:)';
    end
end
Y = double(Y);
height = size(im,1);
width = size(im,2);
clear('im','a');

%-- Split into training and test set
if ~exist('indTs')
   % indTs = [17,22,27,39,45,49]; %[5,11,14];
   indTs = [];
end
indTr = setdiff( 1:size(Y,1),indTs);
Ytr = Y(indTr,:);
Yts = Y(indTs, :);
indPoints = length(indTr);
%---
d = size(Ytr,2);

dataSetName = ['YaleFace' selDir];
%{
    [Y, lbls] = vargplvmLoadData(dataSetName);
    height = lbls(1); width = lbls(2);

    %%
    for i=1:size(Y,1)
        imagesc(reshape(Y(i,:),height, width)), title(num2str(i)),colormap('gray');
        pause
    end
%}

%%



% Set up model
options = vargplvmOptions('dtcvar');
options.kern = mappingKern; %{'rbfard2', 'bias', 'white'};
options.numActive = indPoints;
options.optimiser = 'scg2'; % Set 'scg' if 'scg2' is missing
if ~exist('DgtN') || ~DgtN
    options.enableDgtN = false;
end
options.scaleVal = sqrt(var(Ytr(:)));


% Compute m, the normalised version of Ytr (to be used for
% initialisation of X)
bias = mean(Ytr);
scale = ones(1, d);
if(isfield(options,'scale2var1'))
    if(options.scale2var1)
        scale = std(Ytr);
        scale(find(scale==0)) = 1;
        if(isfield(options, 'scaleVal'))
            warning('Both scale2var1 and scaleVal set for GP');
        end
    end
end
if(isfield(options, 'scaleVal'))
    scale = repmat(options.scaleVal, 1, d);
end
% Remove bias and apply scale.
m = Ytr;
for j = 1:d
    m(:, j) = m(:, j) - bias(j);
    if scale(j)
        m(:, j) = m(:, j)/scale(j);
    end
end

% The length(sharedDims) most important dimensions of the inital X for this
% dataset, are going to be replaced by the shared dimensions of the
% svargplvm model's X.
X_init = ppcaEmbed(m,latentDim);




options.initX = X_init;

% demo using the variational inference method for the gplvm model
fprintf(1,'# Creating the model...\n');

model = vargplvmCreate(latentDim, d, Ytr, options);
% Temporary: in this demo there should always exist the mOrig field
if ~isfield(model, 'mOrig')
    model.mOrig = model.m;
end
model.X = X_init; %%%%%%%
model = vargplvmParamInit(model, model.mOrig, model.X);
model.X = X_init; %%%%%%%

%-- Make sure that the scales for the tied dimensions are at least as large
% as the largest scale.
model.kern.comp{1}.inputScales = invWidthMult./(((max(model.X)-min(model.X))).^2); % Default 5
% For the rbfardjit kernel
model.kern.inputScales = model.kern.comp{1}.inputScales;
params = vargplvmExtractParam(model);
model = vargplvmExpandParam(model, params);
%model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));



model.beta=1/(0.01*var(model.mOrig(:)));
modelInit = model;

capName = dataSetName;
capName(1) = upper(capName(1));
modelType = model.type;
modelType(1) = upper(modelType(1));
fileToSave = ['dem' capName modelType num2str(experimentNo) '.mat'];
%% ---

display = 1;
%     %%%% Optimisation
%     % do not learn beta for few iterations for intitialization
if initVardistIters ~=0
    model.initVardist = 1;
    model.learnBeta = 0; model.learnSigmaf = 0; % This should be merged with the initVardist field
    fprintf(1,'# Intitiliazing the model (fixed beta)...\n');
    model = vargplvmOptimise(model, display, initVardistIters); % Default: 20
    fprintf(1,'1/b = %.4d\n',1/model.beta);
    model.learnSigmaf = 1; model.learnBeta =1; model.initVardist = 0;
end
model.iters = 0;
prunedModelInit = vargplvmPruneModel(modelInit);
clear modelInit

% Optimise the model.
for i=1:length(itNo)
    iters = itNo(i); % default: 2000
    fprintf(1,'\n# Optimising the model for %d iterations (session %d)...\n',iters,i);
    model = vargplvmOptimise(model, display, iters);
    model.iters = model.iters + iters;
    fprintf(1,'1/b = %.4d\n',1/model.beta);
    modelTr = model;
    fprintf(1,'# 1/b=%.4f\n var(m)=%.4f\n',1/model.beta, var(model.mOrig(:)));
    % Save model
    fprintf(1,'# Saving %s\n',fileToSave);
    prunedModel = vargplvmPruneModel(model);
    save(fileToSave, 'prunedModel', 'prunedModelInit');
end
% Just for compatibility.
if strcmp(model.type,'rbfardjit')
    model.kern.comp{1}.inputScales = model.kern.inputScales;
    prunedModel = vargplvmPruneModel(model);
end

save(fileToSave, 'prunedModel', 'prunedModelInit');

model.m = model.mOrig;
%{
sc = 2; % Set to 4 and try again, if you run out of memory
[modelP, newHeight, newWidth] = vargplvmReduceVidModel(model, height, width, sc,sc);
lvmVisualise(modelP, [], 'imageVisualise', 'imageModify', [newHeight newWidth],0,0,1);
clear modelP
figure,bar(model.kern.inputScales)
%}
