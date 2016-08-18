function [Z, varZ, m2] = vargplvmPredictUncertainWrap(model, XextrapMore, YextrapMore, opt)

if nargin < 4,  opt = struct; end
if ~isfield(opt, 'k_it') || isempty(opt.k_it)
    k_it = 100;
else
    k_it = opt.k_it;
end
if ~isfield(opt, 'augmentModel') || isempty(opt.augmentModel)
    opt.augmentModel = true;
end
if ~isfield(opt, 'augmentInd') || isempty(opt.augmentInd)
    opt.augmentInd = false;
end
if ~isfield(opt, 'trCovarsMult') || isempty(opt.trCovarsMult)
    opt.trCovarsMult = 1;
end
if ~isfield(opt, 'ind_it') || isempty(opt.ind_it)
    opt.ind_it = 15;
end
if ~isfield(opt, 'augmentModel_it') || isempty(opt.augmentModel_it)
    opt.augmentModel_it = true;
end




prev = 1;
m2 = model;
Z = zeros(size(YextrapMore));
varZ = zeros(size(YextrapMore));

for jj=1:ceil(size(XextrapMore,1)/k_it)
    fprintf('\n\n# Iteration %d \n', jj)
    XtsCur =  XextrapMore(prev:min(prev+k_it-1, size(XextrapMore,1)),:);
    YCur  =   YextrapMore(prev:min(prev+k_it-1, size(YextrapMore,1)),:);
    [Z(prev:min(prev+k_it-1, size(XextrapMore,1)),:), varZ(prev:min(prev+k_it-1, size(XextrapMore,1)),:)] = ...
        vargplvmPredictUncertainK(m2, XtsCur(1,:), size(XtsCur,1), opt);
    %------ Augment the model with the real data now
    if opt.augmentModel_it
        for i=1:size(XtsCur,1)
            curX = XtsCur(i,:);
            Ypred = YCur(i,:);
            
            m2.N = m2.N + 1;
            %m2.k = m2.k + 1;
            m2.vardist.means = [m2.vardist.means; curX];
            m2.vardist.covars = [m2.vardist.covars; 1e-10*ones(size(curX))];
            
            %m2.vardist.covars = [m2.vardist.covars; curVar+1/m2.beta]; %%%%%%%%%%%%%%%%%%
            
            % Add ind pt if ind_it is true or if it's not logical but specifies
            % every how many steps we need to add the ind. point.
            if (islogical(opt.ind_it) && opt.ind_it)
                m2.X_u = [m2.X_u; curX];
                m2.k = m2.k + size(curX,1);
                pb.printAll = false; %%%%%
            elseif  (~islogical(opt.ind_it) && ~mod(i,opt.ind_it))
                X_u = [m2.X_u; curX];
                res = util_checkCloseMatrixElements(X_u, 0.001);
                if ~ismember(size(X_u,1), res)
                    % Only add ind. pt if it's not close to another one
                    m2.X_u = X_u;
                    m2.k = m2.k + size(curX,1);
                    %pb.nextSymbol = 'o';
                    fprintf('# Adding inducing point...\n');  %%% DEBUG
                else
                    fprintf('# Skipping inducing point...\n');  %%% DEBUG
                end
                pb.printAll = true; %%%%%
            end
            m2.X = [m2.X; curX];
            m2.vardist.numData = m2.N;
            m2.vardist.nParams = 2*prod(size(m2.vardist.means));
            m2.vardist.transforms.index = m2.N*m2.q+1:m2.vardist.nParams;
            m2.numParams = m2.numParams + 2*m2.q;
            m2.nParams = m2.numParams;
            
            % normalize y exactly as model.m is normalized
            my = Ypred - repmat(m2.bias,size(Ypred,1),1);
            my = my./repmat(m2.scale,size(Ypred,1),1);
            
            % change the data (by including the new point and taking only the present indices)
            m2.m = [m2.m; my];
            m2.TrYY = sum(sum(m2.m .* m2.m));
            if isfield(m2, 'fixParamIndices')
                m2.fixParamIndices = 1:2*m2.N*m2.q;
            end
            if isfield(m2, 'fixInducing') && m2.fixInducing
                m2.inducingIndices = 1:m2.k;
            end
            
        end
        %------------
        m2 = vargplvmUpdateStats(m2, m2.X_u);
    end
    prev = prev + k_it;
end
%end
%ZNew3 = Z(1:size(Yextrap,1),:); varZNew3 = varZ(1:size(Yextrap,1),:);