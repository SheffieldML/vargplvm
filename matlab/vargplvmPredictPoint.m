function [X, varX] = vargplvmPredictPoint(dynModel, t_star)

% VARGPLVMPREDICTPOINT Predict the postions of a number of latent points.
% FORMAT
% DESC Given a trained dynamics model dynModel along with a vector of (test) time
% points t_star, this function predicts the latent points corresponding to
% the given time vector. This is done using the conditional GP prior which
% is placed on the latent space and depends on time, i.e. p(X|t). More
% specifically, p(x*|X,t,t*) is just a Gaussian.  So, unlike other
% predictive functions, such as vargplvmOptimisePoint, in this one no
% optimisation is being done because the only given input is time and no
% other information, for example partially observed outputs.
% The number of latent points predicted equals the length of the test time
% vector t_star. 
% ARG dynModel : the dyn. model for which the points will be predicted.
% ARG t_star : the test time points for which the latent points are to
% be predicted. This vector may or may not overlap the time vector for the
% training points.
% RETURN x : the predicted means in the latent space.
% RETURN varx : the predicted variances in the latent space.

% COPYRIGHT: Michalis Titsias, Neil Lawrence, Andreas Damianou 2011
%
% SEEALSO : vargplvmOptimisePoint

% VARGPLVM

N_star = size(t_star,1); % number of test points
K_ss = kernDiagCompute(dynModel.kern, t_star);
K_star = kernCompute(dynModel.kern, dynModel.t, t_star);

X = K_star' * dynModel.vardist.means; % mean

varX = zeros(N_star, dynModel.q); % initialize variances
for q=1:dynModel.q
    invLambda = 1./dynModel.vardist.covars(:,q); 
    Lq = chol(dynModel.Kt + diag(invLambda))'; 
    vq = Lq \ K_star;
    varX(:,q) = K_ss - sum(vq .* vq)';
end
