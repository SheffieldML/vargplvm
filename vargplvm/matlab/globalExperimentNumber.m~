function [expNo, keynames, values] = globalExperimentNumber(demoType, dataType, method, varargin)
% GLOBALEXPERIMENTNUMBER Allows getting/setting global experiment numbers to demos.
% FORMAT
% DESC Allows the user to get/set a global experiment number to each demo.
% This function maintains a file with the experiment numbers used so far
% for each combination of demo and data type.
% ARG method: Alows the user to opt for:
% 'get': get a new number given the demo/data type (in which case an additional argument
% 1 means that the text file is updated to include that, or 0 to not
% include that, default is 1), 
% 'set': specify that the current instance of the demo/data type
% combination will be using the provided experiment number (if it is
% already used for an older experiment, this will be overwritten)
% 'next': equivalent to 'get with an extra argument '0',  meaning that the
% text file is not updated.
% 'free': release as available the provided (in the varargin, as array)
% experiment numbers.

% COPYRIGHT : Andreas C. Damianou, 2011

% VARGPLVM


% constant:
expNumberFile = 'experimentNumbers.txt';

if nargin < 2
    fprintf('Usage: globalExperimentNumber(demoType, dataType, method, varargin)\n');
    return
end
if nargin <3
    method = 'get';
end

[keynames,values]=textread(expNumberFile','%s%s','delimiter','=');

demoType(1) = upper(demoType(1));
dataType(1) = upper(dataType(1));
curKey = ['dem' demoType dataType];

switch method
    case 'set'
        if nargin < 4
            error('With the set method, you must provide the experiment number.');
        end
        expNo = varargin{1};
        indKey = find(strcmp(keynames,curKey));
        if ~isempty(indKey)
            globalExperimentNumber(demoType, dataType, 'free', expNo);
            [keynames,values]=textread(expNumberFile','%s%s','delimiter','=');
        else
            indKey = length(keynames)+1;
            keynames{indKey} = curKey;
            values{indKey} = [];
        end
        expNo = num2str(expNo,'%d ');
        values{indKey} = [values{indKey} ' ' expNo];
        writeBack(expNumberFile, keynames, values);
    case 'get'
        if nargin < 4
            varargin{1} = 1;
        end
        update = varargin{1};
        ind = find(strcmp(keynames,curKey));
        if ~isempty(ind)
            if length(ind)>1
                error('Multiple entries found for key %s!\n',curKey);
            end
            expNo = str2num(values{ind});
            expNo2 = sort(expNo);
            newExp = setdiff(1:expNo2(end), expNo2);
            if isempty(newExp)
                newExp = expNo2(end) + 1;
            else
                newExp = newExp(1);
            end
            expNo = [expNo newExp];
        else
            expNo = 1;
            ind = length(keynames)+1;
        end
        values{ind} = num2str(expNo, '%d ');
        keynames{ind} = curKey;
        if update
             writeBack(expNumberFile, keynames, values);
        end
    case 'next'
        expNo=globalExperimentNumber(demoType, dataType, 'get', 0);
        return
    case 'free'
        %remRange = str2num(varargin{1});
        remRange = varargin{1};
        indKey = find(strcmp(keynames,curKey));
        
        curVal = str2num(values{indKey});
        skipInd = [];
        for i=1:length(remRange)
            skipInd = [skipInd find(curVal == remRange(i))];
        end
        values{indKey} = num2str(setdiff(curVal, curVal(skipInd)),'%d ');
        if isempty(values{indKey})
            inds = setdiff(1:length(keynames),indKey);
            values = values(inds);
            keynames = keynames(inds);
        end
        if ~isempty(skipInd)
            writeBack(expNumberFile, keynames, values);
        end
        expNo = -1;
end


function writeBack(expNumberFile,keynames, values)

fid = fopen(expNumberFile,'w');

for i=1:length(keynames)
    fprintf(fid,[keynames{i} '=' values{i} '\n']);
    fprintf(1,[keynames{i} '=' values{i} '\n']);
end

fclose(fid);
