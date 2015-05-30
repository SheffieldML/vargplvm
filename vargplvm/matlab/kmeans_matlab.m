function [idx, C, sumD, D] = kmeans_matlab(X, k, varargin)
%KMEANS_MATLAB K-means clustering, matlab file copied with a different name
% to avoid conflicts with netlab.
%   IDX = KMEANS(X, K) partitions the points in the N-by-P data matrix
%   X into K clusters.  This partition minimizes the sum, over all
%   clusters, of the within-cluster sums of point-to-cluster-centroid
%   distances.  Rows of X correspond to points, columns correspond to
%   variables.  KMEANS returns an N-by-1 vector IDX containing the
%   cluster indices of each point.  By default, KMEANS uses squared
%   Euclidean distances.
%
%   KMEANS treats NaNs as missing data, and ignores any rows of X that
%   contain NaNs.
%
%   [IDX, C] = KMEANS(X, K) returns the K cluster centroid locations in
%   the K-by-P matrix C.
%
%   [IDX, C, SUMD] = KMEANS(X, K) returns the within-cluster sums of
%   point-to-centroid distances in the 1-by-K vector sumD.
%
%   [IDX, C, SUMD, D] = KMEANS(X, K) returns distances from each point
%   to every centroid in the N-by-K matrix D.
%
%   [ ... ] = KMEANS(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies
%   optional parameter name/value pairs to control the iterative algorithm
%   used by KMEANS.  Parameters are:
%
%   'Distance' - Distance measure, in P-dimensional space, that KMEANS
%      should minimize with respect to.  Choices are:
%          'sqEuclidean'  - Squared Euclidean distance (the default)
%          'cityblock'    - Sum of absolute differences, a.k.a. L1 distance
%          'cosine'       - One minus the cosine of the included angle
%                           between points (treated as vectors)
%          'correlation'  - One minus the sample correlation between points
%                           (treated as sequences of values)
%          'Hamming'      - Percentage of bits that differ (only suitable
%                           for binary data)
%
%   'Start' - Method used to choose initial cluster centroid positions,
%      sometimes known as "seeds".  Choices are:
%          'sample'  - Select K observations from X at random (the default)
%          'uniform' - Select K points uniformly at random from the range
%                      of X.  Not valid for Hamming distance.
%          'cluster' - Perform preliminary clustering phase on random 10%
%                      subsample of X.  This preliminary phase is itself
%                      initialized using 'sample'.
%           matrix   - A K-by-P matrix of starting locations.  In this case,
%                      you can pass in [] for K, and KMEANS infers K from
%                      the first dimension of the matrix.  You can also
%                      supply a 3D array, implying a value for 'Replicates'
%                      from the array's third dimension.
%
%   'Replicates' - Number of times to repeat the clustering, each with a
%      new set of initial centroids.  A positive integer, default is 1.
%
%   'EmptyAction' - Action to take if a cluster loses all of its member
%      observations.  Choices are:
%          'error'     - Treat an empty cluster as an error (the default)
%          'drop'      - Remove any clusters that become empty, and set
%                        the corresponding values in C and D to NaN.
%          'singleton' - Create a new cluster consisting of the one
%                        observation furthest from its centroid.
%
%   'Options' - Options for the iterative algorithm used to minimize the
%       fitting criterion, as created by STATSET.  Choices of STATSET
%       parameters are:
%
%          'Display'  - Level of display output.  Choices are 'off', (the
%                       default), 'iter', and 'final'.
%          'MaxIter'  - Maximum number of iterations allowed.  Default is 100.
%
%   'OnlinePhase' - Flag indicating whether KMEANS should perform an "on-line
%      update" phase in addition to a "batch update" phase.  The on-line phase
%      can be time consuming for large data sets, but guarantees a solution
%      that is a local minimum of the distance criterion, i.e., a partition of
%      the data where moving any single point to a different cluster increases
%      the total sum of distances.  'on' (the default) or 'off'.
%
%   Example:
%
%       X = [randn(20,2)+ones(20,2); randn(20,2)-ones(20,2)];
%       opts = statset('Display','final');
%       [cidx, ctrs] = kmeans(X, 2, 'Distance','city', ...
%                             'Replicates',5, 'Options',opts);
%       plot(X(cidx==1,1),X(cidx==1,2),'r.', ...
%            X(cidx==2,1),X(cidx==2,2),'b.', ctrs(:,1),ctrs(:,2),'kx');
%
%   See also LINKAGE, CLUSTERDATA, SILHOUETTE.

%   KMEANS uses a two-phase iterative algorithm to minimize the sum of
%   point-to-centroid distances, summed over all K clusters.  The first phase
%   uses what the literature often describes as "batch" updates, where each
%   iteration consists of reassigning points to their nearest cluster
%   centroid, all at once, followed by recalculation of cluster centroids.
%   This phase occasionally (especially for small data sets) does not converge
%   to solution that is a local minimum, i.e., a partition of the data where
%   moving any single point to a different cluster increases the total sum of
%   distances.  Thus, the batch phase be thought of as providing a fast but
%   potentially only approximate solution as a starting point for the second
%   phase.  The second phase uses what the literature often describes as
%   "on-line" updates, where points are individually reassigned if doing so
%   will reduce the sum of distances, and cluster centroids are recomputed
%   after each reassignment.  Each iteration during this second phase consists
%   of one pass though all the points.  The on-line phase will converge to a
%   local minimum, although there may be other local minima with lower total
%   sum of distances.  The problem of finding the global minimum can only be
%   solved in general by an exhaustive (or clever, or lucky) choice of
%   starting points, but using several replicates with random starting points
%   typically results in a solution that is a global minimum.
%
% References:
%
%   [1] Seber, G.A.F. (1984) Multivariate Observations, Wiley, New York.
%   [2] Spath, H. (1985) Cluster Dissection and Analysis: Theory, FORTRAN
%       Programs, Examples, translated by J. Goldschmidt, Halsted Press,
%       New York.

%   Copyright 1993-2009 The MathWorks, Inc.
%   $Revision: 1.4.4.11 $  $Date: 2009/05/07 18:31:22 $

if nargin < 2
    error('stats:kmeans:TooFewInputs','At least two input arguments required.');
end

[ignore,wasnan,X] = statremovenan(X);
hadNaNs = any(wasnan);
if hadNaNs
    warning('stats:kmeans:MissingDataRemoved','Ignoring rows of X with missing data.');
end

% n points in p dimensional space
[n, p] = size(X);

pnames = {   'distance'  'start' 'replicates' 'emptyaction' 'onlinephase' 'options' 'maxiter' 'display'};
dflts =  {'sqeuclidean' 'sample'          []         'error'         'on'        []        []        []};
[eid,errmsg,distance,start,reps,emptyact,online,options,maxit,display] ...
    = internal.stats.getargs(pnames, dflts, varargin{:});
if ~isempty(eid)
    error(sprintf('stats:kmeans:%s',eid),errmsg);
end

if ischar(distance)
    distNames = {'sqeuclidean','cityblock','cosine','correlation','hamming'};
    j = strmatch(lower(distance), distNames);
    if length(j) > 1
        error('stats:kmeans:AmbiguousDistance', ...
            'Ambiguous ''Distance'' parameter value:  %s.', distance);
    elseif isempty(j)
        error('stats:kmeans:UnknownDistance', ...
            'Unknown ''Distance'' parameter value:  %s.', distance);
    end
    distance = distNames{j};
    switch distance
        case 'cosine'
            Xnorm = sqrt(sum(X.^2, 2));
            if any(min(Xnorm) <= eps(max(Xnorm)))
                error('stats:kmeans:ZeroDataForCos', ...
                    ['Some points have small relative magnitudes, making them ', ...
                    'effectively zero.\nEither remove those points, or choose a ', ...
                    'distance other than ''cosine''.']);
            end
            X = X ./ Xnorm(:,ones(1,p));
        case 'correlation'
            X = X - repmat(mean(X,2),1,p);
            Xnorm = sqrt(sum(X.^2, 2));
            if any(min(Xnorm) <= eps(max(Xnorm)))
                error('stats:kmeans:ConstantDataForCorr', ...
                    ['Some points have small relative standard deviations, making them ', ...
                    'effectively constant.\nEither remove those points, or choose a ', ...
                    'distance other than ''correlation''.']);
            end
            X = X ./ Xnorm(:,ones(1,p));
        case 'hamming'
            if ~all(ismember(X(:),[0 1]))
                error('stats:kmeans:NonbinaryDataForHamm', ...
                    'Non-binary data cannot be clustered using Hamming distance.');
            end
    end
else
    error('stats:kmeans:InvalidDistance', ...
        'The ''Distance'' parameter value must be a string.');
end

if ischar(start)
    startNames = {'uniform','sample','cluster'};
    j = strmatch(lower(start), startNames);
    if length(j) > 1
        error('stats:kmeans:AmbiguousStart', ...
            'Ambiguous ''Start'' parameter value:  %s.', start);
    elseif isempty(j)
        error('stats:kmeans:UnknownStart', ...
            'Unknown ''Start'' parameter value:  %s.', start);
    elseif isempty(k)
        error('stats:kmeans:MissingK', ...
            'You must specify the number of clusters, K.');
    end
    start = startNames{j};
    if strcmp(start, 'uniform')
        if strcmp(distance, 'hamming')
            error('stats:kmeans:UniformStartForHamm', ...
                'Hamming distance cannot be initialized with uniform random values.');
        end
        Xmins = min(X,[],1);
        Xmaxs = max(X,[],1);
    end
elseif isnumeric(start)
    CC = start;
    start = 'numeric';
    if isempty(k)
        k = size(CC,1);
    elseif k ~= size(CC,1);
        error('stats:kmeans:MisshapedStart', ...
            'The ''Start'' matrix must have K rows.');
    elseif size(CC,2) ~= p
        error('stats:kmeans:MisshapedStart', ...
            'The ''Start'' matrix must have the same number of columns as X.');
    end
    if isempty(reps)
        reps = size(CC,3);
    elseif reps ~= size(CC,3);
        error('stats:kmeans:MisshapedStart', ...
            'The third dimension of the ''Start'' array must match the ''replicates'' parameter value.');
    end
    
    % Need to center explicit starting points for 'correlation'. (Re)normalization
    % for 'cosine'/'correlation' is done at each iteration.
    if isequal(distance, 'correlation')
        CC = CC - repmat(mean(CC,2),[1,p,1]);
    end
else
    error('stats:kmeans:InvalidStart', ...
        'The ''Start'' parameter value must be a string or a numeric matrix or array.');
end

if ischar(emptyact)
    emptyactNames = {'error','drop','singleton'};
    j = strmatch(lower(emptyact), emptyactNames);
    if length(j) > 1
        error('stats:kmeans:AmbiguousEmptyAction', ...
            'Ambiguous ''EmptyAction'' parameter value:  %s.', emptyact);
    elseif isempty(j)
        error('stats:kmeans:UnknownEmptyAction', ...
            'Unknown ''EmptyAction'' parameter value:  %s.', emptyact);
    end
    emptyact = emptyactNames{j};
else
    error('stats:kmeans:InvalidEmptyAction', ...
        'The ''EmptyAction'' parameter value must be a string.');
end

if ischar(online)
    j = strmatch(lower(online), {'on','off'});
    if length(j) > 1
        error('stats:kmeans:AmbiguousOnlinePhase', ...
            'Ambiguous ''OnlinePhase'' parameter value:  %s.', online);
    elseif isempty(j)
        error('stats:kmeans:UnknownOnlinePhase', ...
            'Unknown ''OnlinePhase'' parameter value:  %s.', online);
    end
    online = (j==1);
else
    error('stats:kmeans:InvalidOnlinePhase', ...
        'The ''OnlinePhase'' parameter value must be ''on'' or ''off''.');
end

% 'maxiter' and 'display' are grandfathered as separate param name/value pairs
if ~isempty(display)
    options = statset(options,'Display',display);
end
if ~isempty(maxit)
    options = statset(options,'MaxIter',maxit);
end

options = statset(statset('kmeans'), options);
display = strmatch(lower(options.Display), {'off','notify','final','iter'}) - 1;
maxit = options.MaxIter;

if ~(isscalar(k) && isnumeric(k) && isreal(k) && k > 0 && (round(k)==k))
    error('stats:kmeans:InvalidK', ...
        'X must be a positive integer value.');
    % elseif k == 1
    % this special case works automatically
elseif n < k
    error('stats:kmeans:TooManyClusters', ...
        'X must have more rows than the number of clusters.');
end

% Assume one replicate
if isempty(reps)
    reps = 1;
end

%
% Done with input argument processing, begin clustering
%

dispfmt = '%6d\t%6d\t%8d\t%12g';
if online, Del = NaN(n,k); end % reassignment criterion

totsumDBest = Inf;
emptyErrCnt = 0;
for rep = 1:reps
    switch start
        case 'uniform'
            C = unifrnd(Xmins(ones(k,1),:), Xmaxs(ones(k,1),:));
            % For 'cosine' and 'correlation', these are uniform inside a subset
            % of the unit hypersphere.  Still need to center them for
            % 'correlation'.  (Re)normalization for 'cosine'/'correlation' is
            % done at each iteration.
            if isequal(distance, 'correlation')
                C = C - repmat(mean(C,2),1,p);
            end
            if isa(X,'single')
                C = single(C);
            end
        case 'sample'
            C = X(randsample(n,k),:);
            if ~isfloat(C)      % X may be logical
                C = double(C);
            end
        case 'cluster'
            Xsubset = X(randsample(n,floor(.1*n)),:);
            [dum, C] = kmeans_matlab(Xsubset, k, varargin{:}, 'start','sample', 'replicates',1);
        case 'numeric'
            C = CC(:,:,rep);
    end
    
    % Compute the distance from every point to each cluster centroid and the
    % initial assignment of points to clusters
    D = distfun(X, C, distance, 0);
    [d, idx] = min(D, [], 2);
    m = accumarray(idx,1,[k,1]);
    
    try % catch empty cluster errors and move on to next rep
        
        % Begin phase one:  batch reassignments
        converged = batchUpdate();
        
        % Begin phase two:  single reassignments
        if online
            converged = onlineUpdate();
        end
        
        if ~converged
            warning('stats:kmeans:FailedToConverge', ...
                'Failed to converge in %d iterations%s.',maxit,repsMsg(rep,reps));
        end
        
        % Calculate cluster-wise sums of distances
        nonempties = find(m>0);
        D(:,nonempties) = distfun(X, C(nonempties,:), distance, iter);
        d = D((idx-1)*n + (1:n)');
        sumD = accumarray(idx,d,[k,1]);
        totsumD = sum(sumD);
        
        if display > 1 % 'final' or 'iter'
            disp(sprintf('%d iterations, total sum of distances = %g',iter,totsumD));
        end
        
        % Save the best solution so far
        if totsumD < totsumDBest
            totsumDBest = totsumD;
            idxBest = idx;
            Cbest = C;
            sumDBest = sumD;
            if nargout > 3
                Dbest = D;
            end
        end
        
        % If an empty cluster error occurred in one of multiple replicates, catch
        % it, warn, and move on to next replicate.  Error only when all replicates
        % fail.  Rethrow an other kind of error.
    catch ME
        if reps == 1 || ~isequal(ME.identifier,'stats:kmeans:EmptyCluster')
            rethrow(ME);
        else
            emptyErrCnt = emptyErrCnt + 1;
            warning('stats:kmeans:EmptyCluster', ...
                'Replicate %d terminated: empty cluster created at iteration %d.',rep,iter);
            if emptyErrCnt == reps
                error('stats:kmeans:EmptyClusterAllReps', ...
                    'An empty cluster error occurred in every replicate.');
            end
        end
    end % catch
    
end % replicates

% Return the best solution
idx = idxBest;
C = Cbest;
sumD = sumDBest;
if nargout > 3
    D = Dbest;
end

if hadNaNs
    idx = statinsertnan(wasnan, idx);
end


%------------------------------------------------------------------

    function converged = batchUpdate
        
        % Every point moved, every cluster will need an update
        moved = 1:n;
        changed = 1:k;
        previdx = zeros(n,1);
        prevtotsumD = Inf;
        
        if display > 2 % 'iter'
            disp(sprintf('  iter\t phase\t     num\t         sum'));
        end
        
        %
        % Begin phase one:  batch reassignments
        %
        
        iter = 0;
        converged = false;
        while true
            iter = iter + 1;
            
            % Calculate the new cluster centroids and counts, and update the
            % distance from every point to those new cluster centroids
            [C(changed,:), m(changed)] = gcentroids(X, idx, changed, distance);
            D(:,changed) = distfun(X, C(changed,:), distance, iter);
            
            % Deal with clusters that have just lost all their members
            empties = changed(m(changed) == 0);
            if ~isempty(empties)
                switch emptyact
                    case 'error'
                        error('stats:kmeans:EmptyCluster', ...
                            'Empty cluster created at iteration %d%s.',iter,repsMsg(rep,reps));
                    case 'drop'
                        % Remove the empty cluster from any further processing
                        D(:,empties) = NaN;
                        changed = changed(m(changed) > 0);
                        warning('stats:kmeans:EmptyCluster', ...
                            'Empty cluster created at iteration %d%s.',iter,repsMsg(rep,reps));
                    case 'singleton'
                        warning('stats:kmeans:EmptyCluster', ...
                            'Empty cluster created at iteration %d%s.',iter,repsMsg(rep,reps));
                        
                        for i = empties
                            d = D((idx-1)*n + (1:n)'); % use newly updated distances
                            
                            % Find the point furthest away from its current cluster.
                            % Take that point out of its cluster and use it to create
                            % a new singleton cluster to replace the empty one.
                            [dlarge, lonely] = max(d);
                            from = idx(lonely); % taking from this cluster
                            if m(from) < 2
                                % In the very unusual event that the cluster had only
                                % one member, pick any other non-singleton point.
                                from = find(m>1,1,'first');
                                lonely = find(idx==from,1,'first');
                            end
                            C(i,:) = X(lonely,:);
                            m(i) = 1;
                            idx(lonely) = i;
                            D(:,i) = distfun(X, C(i,:), distance, iter);
                            
                            % Update clusters from which points are taken
                            [C(from,:), m(from)] = gcentroids(X, idx, from, distance);
                            D(:,from) = distfun(X, C(from,:), distance, iter);
                            changed = unique([changed from]);
                        end
                end
            end
            
            % Compute the total sum of distances for the current configuration.
            totsumD = sum(D((idx-1)*n + (1:n)'));
            % Test for a cycle: if objective is not decreased, back out
            % the last step and move on to the single update phase
            if prevtotsumD <= totsumD
                idx = previdx;
                [C(changed,:), m(changed)] = gcentroids(X, idx, changed, distance);
                iter = iter - 1;
                break;
            end
            if display > 2 % 'iter'
                disp(sprintf(dispfmt,iter,1,length(moved),totsumD));
            end
            if iter >= maxit
                break;
            end
            
            % Determine closest cluster for each point and reassign points to clusters
            previdx = idx;
            prevtotsumD = totsumD;
            [d, nidx] = min(D, [], 2);
            
            % Determine which points moved
            moved = find(nidx ~= previdx);
            if ~isempty(moved)
                % Resolve ties in favor of not moving
                moved = moved(D((previdx(moved)-1)*n + moved) > d(moved));
            end
            if isempty(moved)
                converged = true;
                break;
            end
            idx(moved) = nidx(moved);
            
            % Find clusters that gained or lost members
            changed = unique([idx(moved); previdx(moved)])';
            
        end % phase one
        
    end % nested function

%------------------------------------------------------------------

    function converged = onlineUpdate
        
        % Initialize some cluster information prior to phase two
        switch distance
            case 'cityblock'
                Xmid = zeros([k,p,2]);
                for i = 1:k
                    if m(i) > 0
                        % Separate out sorted coords for points in i'th cluster,
                        % and save values above and below median, component-wise
                        Xsorted = sort(X(idx==i,:),1);
                        nn = floor(.5*m(i));
                        if mod(m(i),2) == 0
                            Xmid(i,:,1:2) = Xsorted([nn, nn+1],:)';
                        elseif m(i) > 1
                            Xmid(i,:,1:2) = Xsorted([nn, nn+2],:)';
                        else
                            Xmid(i,:,1:2) = Xsorted([1, 1],:)';
                        end
                    end
                end
            case 'hamming'
                Xsum = zeros(k,p);
                for i = 1:k
                    if m(i) > 0
                        % Sum coords for points in i'th cluster, component-wise
                        Xsum(i,:) = sum(X(idx==i,:), 1);
                    end
                end
        end
        
        %
        % Begin phase two:  single reassignments
        %
        changed = find(m' > 0);
        lastmoved = 0;
        nummoved = 0;
        iter1 = iter;
        converged = false;
        while iter < maxit
            % Calculate distances to each cluster from each point, and the
            % potential change in total sum of errors for adding or removing
            % each point from each cluster.  Clusters that have not changed
            % membership need not be updated.
            %
            % Singleton clusters are a special case for the sum of dists
            % calculation.  Removing their only point is never best, so the
            % reassignment criterion had better guarantee that a singleton
            % point will stay in its own cluster.  Happily, we get
            % Del(i,idx(i)) == 0 automatically for them.
            switch distance
                case 'sqeuclidean'
                    for i = changed
                        mbrs = (idx == i);
                        sgn = 1 - 2*mbrs; % -1 for members, 1 for nonmembers
                        if m(i) == 1
                            sgn(mbrs) = 0; % prevent divide-by-zero for singleton mbrs
                        end
                        Del(:,i) = (m(i) ./ (m(i) + sgn)) .* sum((X - C(repmat(i,n,1),:)).^2, 2);
                    end
                case 'cityblock'
                    for i = changed
                        if mod(m(i),2) == 0 % this will never catch singleton clusters
                            ldist = Xmid(repmat(i,n,1),:,1) - X;
                            rdist = X - Xmid(repmat(i,n,1),:,2);
                            mbrs = (idx == i);
                            sgn = repmat(1-2*mbrs, 1, p); % -1 for members, 1 for nonmembers
                            Del(:,i) = sum(max(0, max(sgn.*rdist, sgn.*ldist)), 2);
                        else
                            Del(:,i) = sum(abs(X - C(repmat(i,n,1),:)), 2);
                        end
                    end
                case {'cosine','correlation'}
                    % The points are normalized, centroids are not, so normalize them
                    normC = sqrt(sum(C.^2, 2));
                    if any(normC < eps(class(normC))) % small relative to unit-length data points
                        error('stats:kmeans:ZeroCentroid', ...
                            'Zero cluster centroid created at iteration %d%s.',iter,repsMsg(rep,reps));
                    end
                    % This can be done without a loop, but the loop saves memory allocations
                    for i = changed
                        XCi = X * C(i,:)';
                        mbrs = (idx == i);
                        sgn = 1 - 2*mbrs; % -1 for members, 1 for nonmembers
                        Del(:,i) = 1 + sgn .*...
                            (m(i).*normC(i) - sqrt((m(i).*normC(i)).^2 + 2.*sgn.*m(i).*XCi + 1));
                    end
                case 'hamming'
                    for i = changed
                        if mod(m(i),2) == 0 % this will never catch singleton clusters
                            % coords with an unequal number of 0s and 1s have a
                            % different contribution than coords with an equal
                            % number
                            unequal01 = find(2*Xsum(i,:) ~= m(i));
                            numequal01 = p - length(unequal01);
                            mbrs = (idx == i);
                            Di = abs(X(:,unequal01) - C(repmat(i,n,1),unequal01));
                            Del(:,i) = (sum(Di, 2) + mbrs*numequal01) / p;
                        else
                            Del(:,i) = sum(abs(X - C(repmat(i,n,1),:)), 2) / p;
                        end
                    end
            end
            
            % Determine best possible move, if any, for each point.  Next we
            % will pick one from those that actually did move.
            previdx = idx;
            prevtotsumD = totsumD;
            [minDel, nidx] = min(Del, [], 2);
            moved = find(previdx ~= nidx);
            if ~isempty(moved)
                % Resolve ties in favor of not moving
                moved = moved(Del((previdx(moved)-1)*n + moved) > minDel(moved));
            end
            if isempty(moved)
                % Count an iteration if phase 2 did nothing at all, or if we're
                % in the middle of a pass through all the points
                if (iter == iter1) || nummoved > 0
                    iter = iter + 1;
                    if display > 2 % 'iter'
                        disp(sprintf(dispfmt,iter,2,nummoved,totsumD));
                    end
                end
                converged = true;
                break;
            end
            
            % Pick the next move in cyclic order
            moved = mod(min(mod(moved - lastmoved - 1, n) + lastmoved), n) + 1;
            
            % If we've gone once through all the points, that's an iteration
            if moved <= lastmoved
                iter = iter + 1;
                if display > 2 % 'iter'
                    disp(sprintf(dispfmt,iter,2,nummoved,totsumD));
                end
                if iter >= maxit, break; end
                nummoved = 0;
            end
            nummoved = nummoved + 1;
            lastmoved = moved;
            
            oidx = idx(moved);
            nidx = nidx(moved);
            totsumD = totsumD + Del(moved,nidx) - Del(moved,oidx);
            
            % Update the cluster index vector, and the old and new cluster
            % counts and centroids
            idx(moved) = nidx;
            m(nidx) = m(nidx) + 1;
            m(oidx) = m(oidx) - 1;
            switch distance
                case 'sqeuclidean'
                    C(nidx,:) = C(nidx,:) + (X(moved,:) - C(nidx,:)) / m(nidx);
                    C(oidx,:) = C(oidx,:) - (X(moved,:) - C(oidx,:)) / m(oidx);
                case 'cityblock'
                    for i = [oidx nidx]
                        % Separate out sorted coords for points in each cluster.
                        % New centroid is the coord median, save values above and
                        % below median.  All done component-wise.
                        Xsorted = sort(X(idx==i,:),1);
                        nn = floor(.5*m(i));
                        if mod(m(i),2) == 0
                            C(i,:) = .5 * (Xsorted(nn,:) + Xsorted(nn+1,:));
                            Xmid(i,:,1:2) = Xsorted([nn, nn+1],:)';
                        else
                            C(i,:) = Xsorted(nn+1,:);
                            if m(i) > 1
                                Xmid(i,:,1:2) = Xsorted([nn, nn+2],:)';
                            else
                                Xmid(i,:,1:2) = Xsorted([1, 1],:)';
                            end
                        end
                    end
                case {'cosine','correlation'}
                    C(nidx,:) = C(nidx,:) + (X(moved,:) - C(nidx,:)) / m(nidx);
                    C(oidx,:) = C(oidx,:) - (X(moved,:) - C(oidx,:)) / m(oidx);
                case 'hamming'
                    % Update summed coords for points in each cluster.  New
                    % centroid is the coord median.  All done component-wise.
                    Xsum(nidx,:) = Xsum(nidx,:) + X(moved,:);
                    Xsum(oidx,:) = Xsum(oidx,:) - X(moved,:);
                    C(nidx,:) = .5*sign(2*Xsum(nidx,:) - m(nidx)) + .5;
                    C(oidx,:) = .5*sign(2*Xsum(oidx,:) - m(oidx)) + .5;
            end
            changed = sort([oidx nidx]);
        end % phase two
        
    end % nested function

end % main function

%------------------------------------------------------------------

function D = distfun(X, C, dist, iter)
%DISTFUN Calculate point to cluster centroid distances.
[n,p] = size(X);
D = zeros(n,size(C,1));
nclusts = size(C,1);

switch dist
    case 'sqeuclidean'
        for i = 1:nclusts
            D(:,i) = (X(:,1) - C(i,1)).^2;
            for j = 2:p
                D(:,i) = D(:,i) + (X(:,j) - C(i,j)).^2;
            end
            % D(:,i) = sum((X - C(repmat(i,n,1),:)).^2, 2);
        end
    case 'cityblock'
        for i = 1:nclusts
            D(:,i) = abs(X(:,1) - C(i,1));
            for j = 2:p
                D(:,i) = D(:,i) + abs(X(:,j) - C(i,j));
            end
            % D(:,i) = sum(abs(X - C(repmat(i,n,1),:)), 2);
        end
    case {'cosine','correlation'}
        % The points are normalized, centroids are not, so normalize them
        normC = sqrt(sum(C.^2, 2));
        if any(normC < eps(class(normC))) % small relative to unit-length data points
            error('stats:kmeans:ZeroCentroid', ...
                'Zero cluster centroid created at iteration %d.',iter);
        end
        
        for i = 1:nclusts
            D(:,i) = max(1 - X * (C(i,:)./normC(i))', 0);
        end
    case 'hamming'
        for i = 1:nclusts
            D(:,i) = abs(X(:,1) - C(i,1));
            for j = 2:p
                D(:,i) = D(:,i) + abs(X(:,j) - C(i,j));
            end
            D(:,i) = D(:,i) / p;
            % D(:,i) = sum(abs(X - C(repmat(i,n,1),:)), 2) / p;
        end
end
end % function

%------------------------------------------------------------------

function [centroids, counts] = gcentroids(X, index, clusts, dist)
%GCENTROIDS Centroids and counts stratified by group.
[n,p] = size(X);
num = length(clusts);
centroids = NaN(num,p);
counts = zeros(num,1);

for i = 1:num
    members = (index == clusts(i));
    if any(members)
        counts(i) = sum(members);
        switch dist
            case 'sqeuclidean'
                centroids(i,:) = sum(X(members,:),1) / counts(i);
            case 'cityblock'
                % Separate out sorted coords for points in i'th cluster,
                % and use to compute a fast median, component-wise
                Xsorted = sort(X(members,:),1);
                nn = floor(.5*counts(i));
                if mod(counts(i),2) == 0
                    centroids(i,:) = .5 * (Xsorted(nn,:) + Xsorted(nn+1,:));
                else
                    centroids(i,:) = Xsorted(nn+1,:);
                end
            case {'cosine','correlation'}
                centroids(i,:) = sum(X(members,:),1) / counts(i); % unnormalized
            case 'hamming'
                % Compute a fast median for binary data, component-wise
                centroids(i,:) = .5*sign(2*sum(X(members,:), 1) - counts(i)) + .5;
        end
    end
end
end % function

%------------------------------------------------------------------

function s = repsMsg(rep,reps)
% Utility for warning and error messages.
if reps == 1
    s = '';
else
    s = sprintf(' during replicate %d',rep);
end
end % function



function [badin,wasnan,varargout]=statremovenan(varargin)
%STATREMOVENAN Remove NaN values from inputs

%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 1.3.2.2 $  $Date: 2005/05/31 16:45:15 $

badin = 0;
wasnan = 0;
n = -1;

% Find NaN, check length, and store outputs temporarily
varargout = cell(nargout,1);
for j=1:nargin
    y = varargin{j};
    if (size(y,1)==1) && (n~=1)
        y =  y';
    end
    
    ny = size(y,1);
    if (n==-1)
        n = ny;
    elseif (n~=ny && ny~=0)
        if (badin==0), badin = j; end
    end
    
    varargout{j} = y;
    
    if (badin==0 && ny>0)
        wasnan = wasnan | any(isnan(y),2);
    end
end

if (badin>0), return; end

% Fix outputs
if (any(wasnan))
    t = ~wasnan;
    for j=1:nargin
        y = varargout{j};
        if (length(y)>0), varargout{j} = y(t,:); end
    end
end
end