function model = vargplvmInitInducing(model, options)
% Sub-sample inducing variables.
model.k = options.numActive;
model.fixInducing = options.fixInducing;


if options.fixInducing
    assert(~isfield(options, 'initX_u') || isempty(options.initX_u));
    assert(model.k == model.N);
    if length(options.fixIndices)~=options.numActive
        %error(['Length of indices for fixed inducing variables must ' ...
        %    'match number of inducing variables']);
        warning(['Length of indices for fixed inducing variables was NOT' ...
            ' equal to the number of inducing variables. This was automatically' ...
            ' fixed by custom fixIndices.']);
        options.fixIndices = 1:options.numActive;
    end
    model.X_u = model.X(options.fixIndices, :);
    model.inducingIndices = options.fixIndices;
else
    if isfield(options, 'initX_u') && ~isempty(options.initX_u)
        fprintf(1,'# Setting inducing points to the given initial value...\n')
        if model.k ~= size(options.initX_u)
            warning('model.k was not equal to the given options.X_u')
            model.k = size(options.initX_u,1);
        end
        model.X_u = options.initX_u;
    else
        if model.k <= model.N
            if ~isfield(options, 'labels')
                if ~isfield(options, 'initX_u')
                    ind = randperm(model.N);
                    ind = ind(1:model.k);
                    model.X_u = model.X(ind, :);
                else
                    model.X_u = options.initX_u;
                end
            else
                % in the case that class labels are supplied, make sure that inducing inputs
                % from all classes are chosen
                [idcs, nSmpls] = class_samples( options.labels, model.k );
                
                count = 1;
                midx = [];
                for inds = idcs
                    ind   = inds{:};
                    ind   = ind(randperm(numel(ind)));
                    idx  = ind(1:nSmpls(count));
                    
                    % test that there is no overlap between index sets
                    assert(isempty(intersect(midx, idx)));
                    midx = [midx, idx];
                    
                    count = count+1;
                end
                model.X_u = model.X(midx,:);
            end
        else
            % Make it work even if k>N
            warning('k > N !')
            % TODO: sample from the variational distr. (this should probably go
            % to the dynamics as well because the vardist. changes in the initialization for the dynamics.
            
            %!!! The following code needs some more testing!
            samplingInd=0; %% TEMP
            if samplingInd
                % This only works if k<= 2*N
                model.X_u=zeros(model.k, model.q);
                ind = randperm(model.N);
                %ind = ind(1:model.N);
                model.X_u(1:model.N,:) = model.X(ind, :);
                
                % The remaining k-N points are sampled from the (k-N) first
                % distributions of the variational distribution (this could be done
                % randomly as well).
                dif=model.k-model.N;
                model.X_u(model.N+1:model.N+dif,:)= model.vardist.means(1:dif,:) + rand(size(model.vardist.means(1:dif,:))).*sqrt(model.vardist.covars(1:dif,:));  % Sampling from a Gaussian.
            else
                model.X_u=zeros(model.k, model.q);
                for i=1:model.k
                    %ind=randi([1 size(model.vardist.means,1)]);
                    % Some versions do not have randi... do it with rendperm
                    % instead:
                    % ceil(size(model.vardist.means,1).*rand) % alternative
                    ind=randperm(size(model.vardist.means,1));
                    ind=ind(1);
                    model.X_u(i,:) = model.vardist.means(ind,:);
                end
            end
        end
    end
    %- Check if some inducing points are too close
    res = util_checkCloseMatrixElements(model.X_u);
    if ~isempty(res)
        warning('The following pairs of inducing points are too close!')
        for ii = 1:size(res,1)%length(res)
            fprintf('%s\n', num2str(res(ii,:)));
        end
    end
end
end

