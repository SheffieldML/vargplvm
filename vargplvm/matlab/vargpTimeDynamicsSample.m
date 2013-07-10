%VARGPTIMEDYNAMICSSAMPLE
function [ySamp, xSamp] = vargpTimeDynamicsSample(model, noiseless)
% VARGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

jitter = 1e-6;


t_star = model.dynamics.t_star;
%-- Precomputations for sampling from q(X_*)
K_ss = kernCompute(model.dynamics.kern, t_star);
K_star = kernCompute(model.dynamics.kern, model.dynamics.t, t_star);


%-- Mean for samples from q(X_*)
X = K_star' * model.dynamics.vardist.means; % mean

xSamp = zeros(size(X));

for q=1:model.dynamics.q
    %-- Covariance for samples from q(X_*)
    invLambda = 1./model.dynamics.vardist.covars(:,q);
    Lq = jitChol(model.dynamics.Kt + diag(invLambda))';
    vq = Lq \ K_star;
    varX = K_ss - vq'*vq; 
        
    %-- Sample from q(X_*) with the mean and covar. calculated above
    try
        xSamp(:,q) = X(:,q) + jitChol(varX)'*randn(size(X(:,q)));
    catch
        xSamp(:,q) = X(:,q) + sqrt(varX)*randn(size(X(:,q)));
    end
end



%-- Sample from p(f_* | X_*)
%- Precomputations:
% Kstar ->  K_{**}
% Kmstar -> K_{M*}

Kstar = kernCompute(model.kern, xSamp);
Kmstar = kernCompute(model.kern, model.X_u, xSamp);

Lm = chol(model.K_uu + jitter*eye(model.k)); %%% in model
invLm = Lm\eye(model.k);   %%% in model
KstarminvL = Kmstar'*invLm;

Ainv = model.P1' * model.P1; % size: NxN
if ~isfield(model,'alpha')
    model.alpha = Ainv*model.Psi1'*model.m; % size: 1xD
end
% mean prediction 
mustar = Kmstar'*model.alpha; % size: 1xD


%- Full covariance
P1Kmstar = model.P1*Kmstar;
CovFstar = Kstar - KstarminvL*KstarminvL' + ...
    (1/model.beta)*(P1Kmstar'*P1Kmstar);

ySamp =  zeros(size(mustar));
D = size(mustar,2);
for d=1:size(mustar,2)
    %- Take the (noise-free) sample
    try
        ySamp(:,d) = mustar(:,d) + jitChol(CovFstar)'*randn(size(mustar,1),1);
    catch
        ySamp(:,d) = mustar(:,d) + sqrt(CovFstar)*randn(size(mustar,1),1);
    end
end

if ~exist('noiseless') || ~noiseless
    %  add also noise
    fprintf(1,'# Adding back the noise...\n')
    ySamp =  ySamp + randn(size(ySamp)).*sqrt(1/model.beta);
end


%%%
% Rescale the mean
ySamp = ySamp.*repmat(model.scale, length(t_star), 1);
% Add the bias back in
ySamp = ySamp + repmat(model.bias, length(t_star), 1);
%%%



% %---------------------
% 
% height = 288; width = 360;
% 
% for i=1:Nstar
%     fr=reshape(ySamp(i,:),height,width);
%     imagesc(fr);
%     colormap('gray');
%     pause(0.2)
% end




