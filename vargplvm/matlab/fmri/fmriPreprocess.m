function [targets, Y] = fmriPreprocess(targetsFile, Y, maskPath)


%--- PART 1: create target's information
%
% Read an fMRI targets file and create a struct 'targets' with all the relevant
% information: targets.runs contains the number of the run for each
% timepoint. The rest of the structure's fields contain the indices of the
% corresponding information.

% open file
if ~exist('targetsFile')
    targetsFile ='chunksTargets_boldDelay3.txt';
end
fid = fopen(targetsFile);

if fid==-1
    error(['Cannot open ' fileName])
end
% parse the file
eof = false;

i=1;
curRun=1;
targets = struct;
targets.runs = [];
targets.binClasses.animals = [];
targets.binClasses.objects = [];
targets.bases = [];
targets.rests = [];
while ~eof
        currentLine = fgetl(fid);
        % check for eof, then CDS. Allow whitespace at the beginning
        if currentLine == -1
            % end of file
            eof = true;
        %elseif ~isempty(regexp(currentLine,'^\s+CDS','match','once'))
        else 
            if strcmp(currentLine(1:2),'a-')
                targets.binClasses.animals = [targets.binClasses.animals i];
                currentLine = currentLine(2:end);
            elseif strcmp(currentLine(1:2),'t-')
                targets.binClasses.objects = [targets.binClasses.objects i];
                currentLine = currentLine(2:end);
            end
            
            if strcmp(currentLine(1:4),'base')
                targets.bases = [targets.bases i];
            elseif strcmp(currentLine(1:4),'rest')
                targets.rests = [targets.rests i];
            end
            targets.runs(i) = str2double(currentLine(end));        
        end
        i = i+1;
end % while eof

% cleanup
fclose(fid);



%--- PART 2: Preprocess the data according to the targets structure.
%
% 

if nargin < 2
    return
end





