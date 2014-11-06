% util_checkCloseMatrixElements Check if some rows of the matrix A are too
% "similar". Similar means their Euclidean distance being smaller than
% tolerance, which if not given is set automatically heuristically
function [res, distances, tolerance] = util_checkCloseMatrixElements(A, tolerance)

C = nan(size(A,1));
for i = 1:size(A,1)
    for j = i+1:size(A,1)
        C(i,j) = dist2(A(i,:), A(j,:));%/norm(abs(A(i,:))+abs(A(j,:)));
    end
end
nn = C(~isnan(C));
distances = C;
C = C./mean(nn);
if nargin < 2 || isempty(tolerance)
    tolerance = 0.0001;
end

%maxDist = dist2(max(A), min(A));
%C = dist2(A, A);
%C = C + eye(size(C))*maxDist*2; % Don't care about diagonal (same point)
%if nargin < 2 || isempty(tolerance)
%    tolerance = (maxDist/mean(sum(A,2).^2))*0.0001;
%    %[~,j]=sort(sum(A,2).^2);
%    %A = A(j,:);
%end
tmp = find(C <= tolerance);
res = [mod(tmp, size(C,1)) ceil(tmp/size(C,1))];
res(res==0) = size(C,1);
res = unique(sort(res,2), 'rows');