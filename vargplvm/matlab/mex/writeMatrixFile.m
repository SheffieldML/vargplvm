function writeMatrixFile(fileName, data)

delete(fileName)
M = [];
sz = size(data);
if length(sz) < 3
    szLimit = size(data,2);
    colName = 'V';
else
    szLimit = size(data,2) * size(data,3);
    colName = 'X';
end

for i=1:szLimit
	if i==1
		M = [M colName num2str(i)];
	else
		M = [M sprintf('\t') colName num2str(i)];
	end
end
dlmwrite(fileName, M, 'delimiter', '')
dlmwrite(fileName, data, '-append','delimiter', '\t')