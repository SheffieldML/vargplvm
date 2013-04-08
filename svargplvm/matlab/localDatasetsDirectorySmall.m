function dir = localDatasetsDirectorySmall

% LOCALDATASETSDIRECTORYSMALL Returns directory where data is stored.
%
%	Description:
%
%	DIR = LOCALDATASETSDIRECTORYSMALL returns the local directory where this file is
%	located, it is assumed that any data sets downloaded have this
%	directory as their base directory. This is for small data files stored
%	locally, e.g. files that can be synced via dropbox etc.
%	 Returns:
%	  DIR - the directory where this m-file is stored.
%	
%
%	See also
%	DATASETSDIRECTORY,DATASETSDIRECTORYLARGE,LVMLOADDATA, MAPLOADDATA
% SHEFFIELDML


% by default return the directory where this file is.
fullSpec = which('localDatasetsDirectorySmall');
ind = max(find(fullSpec == filesep));
dir = fullSpec(1:ind);