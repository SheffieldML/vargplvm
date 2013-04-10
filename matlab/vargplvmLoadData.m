function [Y, lbls, Ytest, lblstest] = vargplvmLoadData(dataset, local, seedVal, field)

% VARGPLVMLOADDATA Load a latent variable model dataset from a local folder
% or from the global repository. This function tries to load the file from
% the following directories (in that order): the local directory where
% the small datasets are kept, the local directory where the large files
% are kept, the global directory as recognised in lvmLoadData, the current
% directory. If the file is found while searching in any of the above
% directories, the function returns so that no more searching is done.
% Searching in local folders succeeds only if such a local folder is added in the
% path and loadLocalData.m is also added in the path (e.g. inside one of
% the local data folders). The loadLocalData.m function has exactly the
% same form as lvmLoadData.m.
% FORMAT
% DESC loads a data set for a latent variable modelling problem.
% ARG dataset : the name of the data set to be loaded.
% ARG seedVal: set a value for the random seeds.
% ARG local: set to true/false for the function to search for dataset in
% the local folders. Default is true.
% RETURN Y : the training data loaded in.
% RETURN lbls : a set of labels for the data (if there are no
% labels it is empty).
% RETURN Ytest : the test data loaded in. If no test set is
% available it is empty.
% RETURN lblstest : a set of labels for the test data (if there are
% no labels it is empty).
%
% Copyright: Andreas C. Damianou, 2011
%
% SEEALSO : lvmLoadData
%

% VARGPLVM

if nargin > 1 && ~isempty(local)
    searchLocally = local;
else
    searchLocally = 1;
end

if nargin > 2 && ~isempty(seedVal)
    randn('seed', seedVal)
    rand('seed', seedVal)
end

if nargin < 3 || isempty(field)
    field = [];
end


% get directory
% if strcmp(fileCategory,'small')
%     baseDir = localDatasetsDirectorySmall;
% elseif strcmp(fileCategory,'large')
%     baseDir = localDatasetsDirectoryLarge;
% end

dirSep = filesep;

% Try the local directories
if searchLocally
    % First, try the local directory (if exists) where the small datasets
    % are kept
    try
        [Y, lbls, Ytest, lblstest]= loadLocalData(dataset, localDatasetsDirectorySmall,dirSep);
        fprintf('# The requested dataset was found in the local directory (for the small files).\n');
        Y = checkField(Y, field);
        return
    catch
        % do nothing
    end

    % If the above fails, try the local directory (if exists) where the
    % large datasets are kept
    try
        [Y, lbls, Ytest, lblstest]=loadLocalData(dataset, localDatasetsDirectoryLarge,dirSep);
        fprintf('# The requested dataset was found in the local directory (for the large files).\n');
        Y = checkField(Y, field);
        return
    catch
        % do nothing
    end
end

% If we reach here nothing of the above contains the dataset. Try the
% global directory.
try
    switch nargout
        case 2
            [Y, lbls] = lvmLoadData(dataset);
        case 3
            [Y, lbls,Ytest] = lvmLoadData(dataset);
        case 4
            [Y, lbls, Ytest, lblstest] = lvmLoadData(dataset);
    end
    fprintf('# The requested dataset was loaded from the global DATASETS directory.\n');
    Y = checkField(Y, field);
    return
catch
    % do nothing
end

% If everything else fails, try to just load the dataset from the current
% dir.
load(dataset);

function Y = checkField(Y, field)
% If Y is a struct, we might only want a specific field to load
if ~isempty(field) && isstruct(Y) && isfield(Y, field)
    Y = Y.(field);
end