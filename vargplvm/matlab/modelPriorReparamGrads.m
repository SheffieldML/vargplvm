function [gVarmeans gVarcovs gDynKern] = modelPriorReparamGrads(dynModel, gVarmeansLik, gVarcovsLik,gInd)
% MODELPRIORREPARAMGRADS Wrapper function for the gradients of the various types of the
% variational GPLVM bound when dynamics is used.
% FORMAT
% DESC provides a wrapper function for the variational GP-LVM, when there
% are dynamics. It takes a dynamics model structure and according to its type it
% calls the appropriate function to calculate the gradient for the
% variational bound. The gradients returned are w.r.t mu_bar and lambda, i.e. the parameters used by the optimiser
% (the ones introduced by the reparametrization), not the original ones (mu and S).
% These gradients are obtained in two steps: firstly, another funtion is used to calculate the
% gradients that correspond only to the likelihood part of the bound and
% are w.r.t the original parameters mu and S. These quantities are the arguments
% gVarmenasLik and gVarcovsLik and are calculated using the already
% implemented code for the static vargplvm, which assumes that the original
% parameters are not coupled.
% The current funtions receives these quantities and a) amends
% with partial derivatives because the variational parameters are coupled
% via Kt for the dyn. gplvm b) Computes the whole gradients for both parts
% of the bound, the one corresponding to the likelihood and the one
% corresponding to the prior. This must be done in a single function
% because the final formula contains both parts in a nonlinear form.
% c) Also the kernel hyperparameters for the dynamics kernel are being retured.
%
% See the dyn. vargplvm notes for more details.
%
% ARG dynModel : the dynamics model structure for which the gradients are
% to be computed.
% ARG gVarmeansLik, gVarcovsLik: the gradients for the VAR-GPLVM model computed
% for the original parameters and only for the likelihood term.
% ARG gInd: in case the inducing points are tied to the variational means
% this is the partial derivatives of the likelihood part of the variational
% bound w.r.t the inducing points X_u, otherwise it is just [].
% RETURN gVarmeans, gVarcovs : the gradients for the "reparametrized" means
% and covariances (the ones that are visible to the optimiser, not the
% original ones) for the VAR-GPLVM model.
% RETURN gDynKern: the gradient w.r.t the hyperparameters of the dynamics
% kernel
% 
% SEEALSO : vargplvmLogLikeGradients vargpTimeDynamicsPriorReparamGrads
%
% COPYRIGHT : Michalis K. Titsias, 2010-2011
% COPYRIGHT : Neil D. Lawrence, 2010-2011
% COPYRIGHT : Andreas C. Damianou, 2010-2011

% VARGPLVM


fhandle = str2func([dynModel.type 'PriorReparamGrads']);
[gVarmeans gVarcovs gDynKern] = fhandle(dynModel, gVarmeansLik, gVarcovsLik, gInd);


